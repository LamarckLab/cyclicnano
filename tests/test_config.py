"""Config loading, derived values and startup validation."""
from __future__ import annotations

from pathlib import Path

from cyclicnano.config import load_config, parse_overrides
from cyclicnano.filters import FilterSet
from cyclicnano.geometry import backbone_metrics
from synthetic import cyclic_oligomer


def _cfg(*overrides):
    return load_config(overrides=parse_overrides(list(overrides)))


def test_symmetry_order_parsed():
    assert _cfg("target.symmetry=C5").sym_order == 5
    assert _cfg("target.symmetry=c3").sym_order == 3


def test_unsupported_symmetry_rejected():
    try:
        _cfg("target.symmetry=D2").sym_order
        raise AssertionError("expected a ValueError")
    except ValueError:
        pass


def test_contig_length_is_total_not_per_chain():
    """RFdiffusion's symmetric config takes the total across all chains."""
    assert _cfg("target.symmetry=C5", "target.monomer_length=60").total_length == 300
    assert _cfg("target.symmetry=C3", "target.monomer_length=80").total_length == 240


def test_overrides_are_nested_and_typed():
    cfg = _cfg("design.n_seq_per_backbone=32", "compute.gpus=[2,3]")
    assert cfg.get("design.n_seq_per_backbone") == 32
    assert cfg.get("compute.gpus") == [2, 3]


def test_malformed_override_rejected():
    try:
        parse_overrides(["target.symmetry"])
        raise AssertionError("expected a ValueError")
    except ValueError:
        pass


def test_default_config_validates_cleanly():
    assert _cfg("target.symmetry=C5", "target.monomer_length=60").validate() == []


def test_untied_chains_rejected():
    problems = _cfg("design.tie_chains=false").validate()
    assert any("tie_chains" in p for p in problems)


def test_msa_mode_other_than_single_sequence_is_flagged():
    problems = _cfg("fold.msa_mode=mmseqs2_uniref_env").validate()
    assert any("msa_mode" in p for p in problems)


def test_null_gpus_is_allowed_and_leaves_the_card_choice_alone():
    from cyclicnano.runner import Runner
    assert _cfg("compute.gpus=null").validate() == []
    assert "CUDA_VISIBLE_DEVICES" not in Runner(kind="local", gpus=[]).env_vars()


def test_explicit_gpu_list_still_pins():
    from cyclicnano.runner import Runner
    env = Runner(kind="local", gpus=[2, 3]).env_vars()
    assert env["CUDA_VISIBLE_DEVICES"] == "2,3"
    assert env["CUDA_DEVICE_ORDER"] == "PCI_BUS_ID"


def test_implausible_monomer_length_rejected():
    assert any("monomer_length" in p for p in _cfg("target.monomer_length=5").validate())


def test_require_raises_for_missing_key():
    try:
        _cfg().require("no.such.key")
        raise AssertionError("expected a KeyError")
    except KeyError:
        pass


def test_shipped_default_rules_parse_against_real_metrics():
    cfg = _cfg("target.symmetry=C5")
    metrics = backbone_metrics(cyclic_oligomer(n_sym=5), expected_sym=5)
    FilterSet("backbone", cfg.get("backbone_filter.rules")).check_syntax(metrics)


def test_shipped_validate_rules_parse():
    cfg = _cfg()
    sample = {"success_rate": 0.5, "best_rmsd": 0.8, "median_rmsd": 1.1,
              "n_seqs": 10, "n_seqs_folded": 10, "best_plddt": 90.0}
    FilterSet("validate", cfg.get("validate.rules")).check_syntax(sample)


def test_zero_seed_is_rejected():
    """ProteinMPNN reads --seed 0 as 'pick a random one', which breaks reproducibility."""
    assert any("run.seed" in p for p in _cfg("run.seed=0").validate())
    assert _cfg("run.seed=1").validate() == []


# ---------------------------------------------------------------- layered run configs
def test_several_config_files_merge_in_order():
    """Later files win, so a variant file need only state what differs from the base."""
    cfg = load_config(["configs/runs/base.yaml", "configs/runs/c3.yaml"])
    assert cfg.symmetry == "C3"                      # from the variant
    assert cfg.get("target.n_backbones") == 500      # from the base
    assert cfg.get("design.n_seq_per_backbone") == 20


def test_shipped_run_configs_differ_only_in_symmetry_and_naming():
    """The point of the split is that everything else cannot drift apart."""
    import yaml
    from cyclicnano.config import REPO_ROOT
    c5 = yaml.safe_load((REPO_ROOT / "configs/runs/c5.yaml").read_text(encoding="utf-8"))
    c3 = yaml.safe_load((REPO_ROOT / "configs/runs/c3.yaml").read_text(encoding="utf-8"))
    assert set(c5) == set(c3) == {"run", "target"}
    assert set(c5["target"]) == {"symmetry"}
    assert c5["target"]["symmetry"] == "C5" and c3["target"]["symmetry"] == "C3"


def test_both_shipped_run_configs_validate():
    for variant in ("c5", "c3"):
        cfg = load_config(["configs/runs/base.yaml", f"configs/runs/{variant}.yaml"])
        assert cfg.validate() == [], (variant, cfg.validate())


def test_contig_length_follows_the_symmetry_of_each_variant():
    assert load_config(["configs/runs/base.yaml", "configs/runs/c5.yaml"]).total_length == 300
    assert load_config(["configs/runs/base.yaml", "configs/runs/c3.yaml"]).total_length == 180


def test_contig_total_is_divisible_by_the_symmetry_order():
    """RFdiffusion rejects a symmetric contig that does not divide evenly."""
    for variant, order in (("c5", 5), ("c3", 3)):
        cfg = load_config(["configs/runs/base.yaml", f"configs/runs/{variant}.yaml"])
        assert cfg.total_length % order == 0, (variant, cfg.total_length)


def test_a_single_config_path_still_works():
    """The list form must not break callers that pass one file."""
    assert load_config("configs/runs/c3.yaml").symmetry == "C3"


# ---------------------------------------------------------------- provenance
def test_resolved_yaml_records_every_source():
    cfg = load_config(["configs/runs/base.yaml", "configs/runs/c5.yaml"],
                      "configs/profiles/amax.yaml",
                      parse_overrides(["target.n_backbones=7"]))
    text = cfg.resolved_yaml()
    for fragment in ("default.yaml", "base.yaml", "c5.yaml", "amax.yaml", "command line"):
        assert fragment in text, fragment


def test_resolved_yaml_reloads_to_the_same_config():
    """The dump must be a usable config, not just a readable one."""
    import tempfile
    cfg = load_config(["configs/runs/base.yaml", "configs/runs/c3.yaml"],
                      overrides=parse_overrides(["target.n_backbones=7"]))
    tmp = Path(tempfile.mkdtemp(prefix="cyclicnano_resolved_")) / "resolved.yaml"
    tmp.write_text(cfg.resolved_yaml(), encoding="utf-8")
    again = load_config(str(tmp))
    assert again.symmetry == "C3"
    assert again.get("target.n_backbones") == 7      # the override survived the round trip
    assert again.data == cfg.data


# ---------------------------------------------------------------- output location
def test_outdir_is_outroot_plus_name():
    cfg = load_config(overrides=parse_overrides(
        ["run.outroot=/data/runs", "run.name=c5", "run.outdir=null"]))
    assert cfg.outdir.as_posix().endswith("/data/runs/c5")


def test_explicit_outdir_overrides_outroot():
    cfg = load_config(overrides=parse_overrides(
        ["run.outroot=/data/runs", "run.name=c5", "run.outdir=/elsewhere/here"]))
    assert cfg.outdir.as_posix().endswith("/elsewhere/here")


def test_posix_outroot_is_absolute_even_on_windows():
    """A config written for a Linux host must not be resolved against the repo root."""
    from cyclicnano.config import REPO_ROOT
    cfg = load_config(overrides=parse_overrides(
        ["run.outroot=/data/lmk/cyclicnano_outputs", "run.name=c3", "run.outdir=null"]))
    # is_relative_to, not a substring test: after the repository was renamed to
    # cyclicnano its path became a prefix of cyclicnano_outputs, and a naive
    # "not in" check started failing on a layout that is perfectly correct.
    assert not cfg.outdir.is_relative_to(REPO_ROOT)
    assert cfg.outdir.as_posix() == "/data/lmk/cyclicnano_outputs/c3"


def test_relative_outroot_resolves_against_the_repository():
    from cyclicnano.config import REPO_ROOT
    cfg = load_config(overrides=parse_overrides(
        ["run.outroot=runs", "run.name=c5", "run.outdir=null"]))
    assert cfg.outdir == (REPO_ROOT / "runs" / "c5").resolve()


def test_amax_profile_sends_results_outside_the_code_directory():
    from cyclicnano.config import REPO_ROOT
    for variant in ("c5", "c3"):
        cfg = load_config(["configs/runs/base.yaml", f"configs/runs/{variant}.yaml"],
                          "configs/profiles/amax.yaml")
        assert cfg.outdir.as_posix() == f"/data/lmk/cyclicnano_outputs/{variant}"
        assert not cfg.outdir.is_relative_to(REPO_ROOT)


def test_shipped_configs_define_no_key_default_yaml_does_not():
    """The regression this exists for: an edit dropped the `generate:` header from
    base.yaml, so batch_size and noise_scale sat under `target:` where nothing reads
    them. Every test passed, the config said noise 0.5, and the run used 1.0."""
    for variant in ("c5", "c3"):
        cfg = load_config(["configs/runs/base.yaml", f"configs/runs/{variant}.yaml",
                           "configs/presets/pure_three_helix.yaml"],
                          "configs/profiles/amax.yaml")
        assert cfg.unknown_keys() == [], (variant, cfg.unknown_keys())
        assert cfg.validate() == [], (variant, cfg.validate())


def test_a_misplaced_key_is_reported_rather_than_ignored():
    cfg = load_config(overrides=parse_overrides(["target.noise_scale=0.5"]))
    assert "target.noise_scale" in cfg.unknown_keys()
    assert any("target.noise_scale" in p for p in cfg.validate())


def test_free_form_subtrees_are_not_policed():
    """backend.tool_env is an empty dict in the defaults by design: what goes in it is
    machine-specific, so the check must stop descending there rather than reject it."""
    cfg = load_config(profile="configs/profiles/amax.yaml")
    assert cfg.get("backend.tool_env.colabfold")          # the profile does fill it
    assert cfg.unknown_keys() == []


def test_the_shipped_base_actually_delivers_the_noise_scale():
    """Reads the value through the same path the stage does, then checks it reaches
    the emitted command. Asserting on the YAML alone would not have caught the bug."""
    from cyclicnano.s00_generate import build_command
    cfg = load_config(["configs/runs/base.yaml", "configs/runs/c5.yaml"])
    assert cfg.get("generate.noise_scale") == 0.5
    argv = [str(a) for a in build_command(cfg, Path("/out/design"), 10)]
    assert "denoiser.noise_scale_ca=0.5" in argv
    assert "denoiser.noise_scale_frame=0.5" in argv


def test_shipped_configs_are_present_and_tracked():
    """They were not tracked: an unanchored runs/ in .gitignore also matched
    configs/runs/, so git add skipped them without a word.

    Nothing failed locally, because the files existed in the working tree. CI
    checked out a tree without them and every config test broke.
    """
    import subprocess
    from cyclicnano.config import REPO_ROOT
    required = ["configs/default.yaml", "configs/runs/base.yaml",
                "configs/runs/c5.yaml", "configs/runs/c3.yaml",
                "configs/profiles/amax.yaml"]

    for name in required:
        assert (REPO_ROOT / name).exists(), f"{name} is missing from the tree"

    if not (REPO_ROOT / ".git").exists():
        return                                   # an exported tree has nothing to check
    tracked = subprocess.run(["git", "ls-files", "configs/"], cwd=REPO_ROOT,
                             capture_output=True, text=True).stdout.split()
    for name in required:
        assert name in tracked, f"{name} exists but is not under version control"


def test_filter_preset_layers_onto_a_run_config():
    """The filter is applied by stacking a third file, not by editing the run config."""
    cfg = load_config(["configs/runs/base.yaml", "configs/runs/c5.yaml",
                       "configs/presets/pure_three_helix.yaml"])
    assert cfg.symmetry == "C5"                              # variant survives
    assert cfg.get("target.n_backbones") == 500              # base survives
    assert cfg.get("backbone_filter.rules") == [             # preset wins on its own key
        "strand_frac == 0", "n_helices == 3"]


def test_every_shipped_preset_parses_against_the_metric_set():
    """A preset that names a metric the pipeline does not compute fails at startup,
    which is far too late if the run was launched overnight."""
    from cyclicnano.cli import _SAMPLE_METRICS
    from cyclicnano.config import REPO_ROOT
    from cyclicnano.filters import FilterSet

    presets = sorted((REPO_ROOT / "configs" / "presets").glob("*.yaml"))
    assert presets, "no presets shipped"
    for preset in presets:
        cfg = load_config(["configs/runs/base.yaml", "configs/runs/c5.yaml",
                           f"configs/presets/{preset.name}"])
        rules = cfg.get("backbone_filter.rules")
        assert rules, f"{preset.name} sets no rules"
        FilterSet("backbone", rules).check_syntax(_SAMPLE_METRICS)
