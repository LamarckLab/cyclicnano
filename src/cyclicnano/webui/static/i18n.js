/* The two languages the page speaks.
 *
 * English is the default. Chinese is a translation of the prose only: every
 * technical term keeps its usual spelling, because RMSD, pLDDT, ddG, REU,
 * RFdiffusion, ProteinMPNN, AlphaFold2, Rosetta, FASTA and PDB are what they are
 * called in the lab and in the literature, and translating them would make the
 * page harder to read rather than easier.
 *
 * `{0}`, `{1}` are positional. A key missing from `zh` falls back to `en`, so a
 * half-finished translation shows English rather than a key name.
 */
'use strict';

const STRINGS = {
  en: {
    /* --------------------------------------------------------- the projects */
    'projects.title': 'Projects',
    'projects.sub': 'One project is one run directory: the design target, the stages '
      + 'that have run on it, and everything they produced. Pick one to see the '
      + 'backbones it delivered.',
    'projects.new': 'New project',
    'projects.none': 'No projects yet',
    'projects.none.detail': 'Nothing under the output root looks like a project. A '
      + 'directory counts as one once it holds a project.json, a tables/backbones.csv '
      + 'or an 04_delivery/delivery.json. Directories whose name starts with an '
      + 'underscore are never listed.',
    'projects.empty': 'nothing computed yet',
    'projects.empty.hint': 'This project has a name and a target and no results. '
      + 'Submit a job to start the pipeline on it.',
    'projects.mismatch': 'declared {0} but the tables measured {1}',
    'projects.created': 'created {0}',

    'new.title': 'New project',
    'new.sub': 'A project is a directory under the output root. Creating one writes '
      + 'its name and target and nothing else: naming a target and spending GPU hours '
      + 'on it are separate decisions, so no stage runs here.',
    'new.name': 'Directory name',
    'new.name.help': 'Becomes <output root>/<name>/ and part of the page URL. Letters, '
      + 'digits, underscore, dot and hyphen.',
    'new.label': 'Title',
    'new.label.help': 'Shown on the card. Free text; defaults to the directory name.',
    'new.symmetry': 'Symmetry',
    'new.symmetry.help': 'The n of Cn. Recorded with the project and used as the '
      + 'default when a job is submitted.',
    'new.monomer': 'Monomer length (aa)',
    'new.monomer.help': 'Residues per chain. The RFdiffusion contig total becomes n '
      + 'times this.',
    'new.note': 'Note',
    'new.note.help': 'Optional. What this project is for, in your own words.',
    'new.create': 'Create the project',
    'new.cancel': 'Back',
    'new.disabled': 'Project creation is not enabled on this service.',

    /* ------------------------------------------------- the job form's fields */
    'jobs.field.name': 'Run name',
    'jobs.field.name.help': 'Results go to <output root>/<name>/. An existing name '
      + 'resumes that run, skipping the stages already done.',
    'jobs.field.symmetry': 'Symmetry',
    'jobs.field.symmetry.help': 'The n of Cn. The contig total becomes n times the '
      + 'monomer length.',
    'jobs.field.monomer_length': 'Monomer length (aa)',
    'jobs.field.n_backbones': 'Backbones to generate',
    'jobs.field.n_backbones.help': 'The stage-01 pass rate is about 2.5% for C5 and 9% '
      + 'for C3, so this number has to be large.',
    'jobs.field.noise_scale': 'RFdiffusion noise',
    'jobs.field.noise_scale.help': '0.5 yields twice as many three-helix backbones as '
      + '1.0. 0.0 yields more helix but mostly long hairpins, and was ruled out on '
      + 'measurement.',
    'jobs.field.batch_size': 'Backbones per batch',
    'jobs.field.batch_size.help': 'Batching is only for resuming; the batch boundaries '
      + 'are not themselves reproducible.',
    'jobs.field.n_seq_per_backbone': 'Sequences per backbone',
    'jobs.field.rmsd_cutoff': 'RMSD cutoff (\u00c5)',
    'jobs.field.rmsd_cutoff.help': 'Whether a sequence folded back into the designed '
      + 'subunit. The literature usually says 2.0; this pipeline was asked for 1.0.',
    'jobs.field.plddt_cutoff': 'pLDDT cutoff',
    'jobs.field.seed': 'Random seed',
    'jobs.submit.here': 'Submit a job',
    'jobs.submit.locked': 'The run name and symmetry come from this project.',

    /* ------------------------------------------------- the interface residues */
    'iface.title': 'Interface residues',
    'iface.sub': '{0} residues, heavy atoms within {1} \u00c5 across the two chains',
    'iface.chain': 'chain {0}',
    'iface.none': 'No contacts found between the two chains.',
    'iface.loading': 'Finding the interface\u2026',
    'iface.fail': 'Could not read the interface: {0}',

    /* ----------------------------------------------------------- the chrome */
    'nav.submit': 'Submit',
    'theme.toggle': 'Switch between light and dark',
    'lang.toggle': 'Switch language',
    'lang.name': '中文',
    'root.unreadable': ' (unreadable)',

    /* ------------------------------------------------------- level 1: runs */
    'runs.complete': 'complete',
    'runs.incomplete': 'not marked complete',
    'runs.undelivered': 'no delivery yet',
    'runs.undelivered.title': 'Has tables/ but no 04_delivery/',
    'runs.rate': 'overall yield {0}%',
    'runs.interface': '{0} sequences scored for interface',
    'runs.updated': 'updated {0}',

    'funnel.generated': 'generated',
    'funnel.generated.title': 'Symmetric backbones RFdiffusion produced',
    'funnel.filtered': 'passed filter',
    'funnel.filtered.title': 'Backbones that passed the stage-01 shape filter',
    'funnel.sequences': 'MPNN sequences',
    'funnel.sequences.title': 'Sequences ProteinMPNN designed for the backbones that passed',
    'funnel.predictions': 'AF2 predictions',
    'funnel.predictions.title': 'Sequences whose no-MSA monomer prediction completed',
    'funnel.folded': 'folded back',
    'funnel.folded.title': 'Sequences that met both the RMSD and the pLDDT cutoff',
    'funnel.delivered': 'delivered',
    'funnel.delivered.title': 'Backbones whose success rate met the cutoff',
    'funnel.kept': ', {0}% of the step before',

    /* --------------------------------------------------- level 2: one run */
    'run.sub': '{0} backbones passed the filter, {1} of them met the success-rate '
      + 'cutoff{2}{3}. Click a row for its sequences and structures; the column '
      + 'headers sort.',
    'run.sub.monomer': ' · {0} aa monomer',
    'run.nodelivery.title': '{0} · no delivery table yet',
    'run.nodelivery.sub': 'This run has not reached 04_delivery, so there is no '
      + 'per-backbone table. What each stage dropped is in the funnel on the runs page.',
    'run.nodelivery.hint': 'To fill the table in, submit the same run name again with '
      + 'the stage range set to 04 or 05. The stages that are done are skipped, so no '
      + 'GPU work is repeated.',

    'col.rank': 'Rank',
    'col.rank.title': 'Delivery rank; a backbone below the cutoff has none',
    'col.backbone': 'Backbone',
    'col.delivered': 'Delivered',
    'col.delivered.title': 'Whether the monomer success rate met the cutoff',
    'col.success': 'Success rate',
    'col.success.title': 'Share of this backbone’s sequences that met both the '
      + 'RMSD and the pLDDT cutoff',
    'col.rmsd.dist.title': 'Every sequence of this backbone by CA RMSD to the design, '
      + 'as a density curve over a shared 0–4 Å axis with the real values as '
      + 'ticks beneath it. The numbers at the ends are this backbone’s own minimum '
      + 'and maximum; the dividing line is the cutoff. Height is scaled per row, so a '
      + 'taller peak means more concentrated here, not more sequences.',
    'col.plddt.dist.title': 'Every sequence of this backbone by mean pLDDT, as a '
      + 'density curve over a shared 50–100 axis with the real values as ticks '
      + 'beneath it. The numbers at the ends are this backbone’s own minimum and '
      + 'maximum; the dividing line is the cutoff.',
    'col.ddg.dist.title': 'The sequences that were scored, by interface energy, as a '
      + 'density curve over a shared −85–+5 REU axis with the real values '
      + 'as ticks beneath it. The numbers at the ends are this backbone’s own '
      + 'minimum and maximum; the dividing line is the cutoff. Sequences that were '
      + 'never scored are left out and counted in the tooltip.',
    'col.iface.rate': 'Interface rate',
    'col.iface.rate.title': 'Share of the scored sequences below the ddG cutoff. '
      + 'Reported beside the monomer success rate, never folded into it.',
    'col.seqs': 'Seqs',
    'col.seqs.title': 'Sequences ProteinMPNN designed',
    'col.preds': 'Preds',
    'col.preds.title': 'Sequences whose AF2 monomer prediction completed',
    'col.failreason': 'Fail reason',
    'col.index': '#',
    'col.pass': 'Pass',
    'col.pass.title': 'Whether both the RMSD and the pLDDT cutoff were met',
    'col.rmsd': 'RMSD',
    'col.rmsd.title': 'CA RMSD between the predicted monomer and the design, after '
      + 'superposition (Å)',
    'col.plddt': 'pLDDT',
    'col.plddt.title': 'Mean pLDDT of the predicted monomer',
    'col.ddg': 'ddG',
    'col.ddg.title': 'Labelled ddG because that is what the group calls it. The '
      + 'number is Rosetta InterfaceAnalyzerMover’s dG_separated: the score of '
      + 'the complex minus the two separated and repacked chains, in REU — not '
      + 'kcal/mol. Strictly it is a binding energy, not a ΔΔG, which is '
      + 'the difference between two of them and needs a reference state that a de '
      + 'novo design does not have. It decides interface_pass and nothing else.',
    'col.ddg.open': 'Click to show the relaxed dimer this energy was measured on, with its interface residues',
    'col.iface': 'Interface',
    'col.iface.title': 'Whether this sequence is below the ddG cutoff. Describes the '
      + 'interface; filters nothing.',
    'col.sequence': 'Sequence',

    /* --------------------------------------------------------- the export */
    'export.scope': 'Export',
    'export.fasta': 'Download FASTA',
    'export.csv': 'Download CSV (opens in Excel)',
    'export.scope.delivered': 'Delivered backbones · passing sequences',
    'export.scope.delivered.help': 'Passing sequences from backbones whose success '
      + 'rate met the cutoff. This is the set you would order synthesis from.',
    'export.scope.passed': 'All passing sequences',
    'export.scope.passed.help': 'Plus the sequences that folded back on a backbone '
      + 'whose overall rate did not — where a lucky single design hides.',
    'export.scope.all': 'Every sequence, failures included',
    'export.scope.all.help': 'Everything in the delivery, for the record.',
    'export.note': 'The energy column in the exported files is named dG_interface; it '
      + 'is the same number the table shows as {0}. The files use the stored field '
      + 'name so they stay traceable once they leave this page.',

    /* ------------------------------------------------- level 3: a backbone */
    'detail.sequences': 'Sequences',
    'detail.structure': 'Structure',
    'detail.criteria': 'A sequence passes at RMSD ≤ {0} Å and mean pLDDT '
      + '≥ {1}; those two decide the success rate and the delivery verdict. The '
      + 'interface is judged separately, at {2} < {3} REU, and gates nothing.',
    'detail.criteria.noddg': 'A sequence passes at RMSD ≤ {0} Å and mean '
      + 'pLDDT ≥ {1}. {2} describes the interface and gates nothing.',
    'fact.success': 'Success rate',
    'fact.passed': 'Passed',
    'fact.bestrmsd': 'Best RMSD',
    'fact.bestplddt': 'Best pLDDT',
    'fact.iface': 'Interface rate',
    'fact.ifacecount': '{0} of {1} scored',
    'fact.delivered': 'Delivered',
    'yes': 'yes',
    'no': 'no',
    'detail.pickrow': 'Click a row to show that sequence’s structure.',
    'detail.best': 'best',

    /* --------------------------------------------------------- the viewer */
    'viewer.structure': 'Structure',
    'viewer.colour': 'Colour',
    'viewer.group.backbone': 'Design backbone',
    'viewer.group.current': 'Current sequence {0}',
    'viewer.group.other': 'Other sequences',
    'viewer.missing': ' (file missing)',
    'struct.backbone': 'design backbone',
    'struct.monomer': 'AF2 monomer',
    'struct.oligomer': '{0}-mer, design geometry',
    'struct.oligomer.3': 'trimer, design geometry',
    'struct.oligomer.5': 'pentamer, design geometry',
    'struct.dimer': 'relaxed dimer',
    'note.backbone': 'The filtered RFdiffusion backbone: main chain only, no side chains.',
    'note.monomer': 'The no-MSA monomer prediction. The B-factor column holds pLDDT. '
      + 'AlphaFold2 writes its own coordinate frame, so this does not sit on the '
      + 'design backbone and the RMSD is measured after superposition.',
    'note.oligomer': 'The monomer prediction placed on every subunit of the design '
      + 'backbone, assembled here for viewing. Rigid-body superposition with no '
      + 'relaxation, so neighbouring subunits may clash. Nothing in this pipeline '
      + 'predicts the assembly: this is the geometry the design asks for, carrying '
      + 'the side chains the sequence implies, and not evidence that the sequence '
      + 'assembles into it.',
    'note.dimer': 'Two copies of the monomer prediction superposed onto chains A and B '
      + 'of the design backbone, then relaxed under restraints. {0} is measured here, '
      + 'so it describes the interface that was designed. The relax displaces the whole '
      + 'object by 2–3 Å, but changes its internal geometry by only 0.8 Å '
      + 'and the distance between the two subunits by under 1 Å: the drift is rigid '
      + 'and the interface geometry held. The residues listed below are marked in '
      + 'the viewer as selected; click empty space in the viewport to clear them.',
    'colour.chain': 'by chain',
    'colour.plddt': 'by pLDDT',
    'colour.ss': 'by secondary structure',
    'colour.na': ' (not available here)',
    'legend.plddt': 'pLDDT',
    'legend.plddt.vhigh': '> 90 very high',
    'legend.plddt.high': '70–90 confident',
    'legend.plddt.low': '50–70 low',
    'legend.plddt.vlow': '≤ 50 very low',
    'legend.plddt.source': '(AlphaFold’s own bands, read from the B-factor column)',
    'legend.molstar': 'Colours come from Mol*’s {0} theme.',
    'viewer.download': 'Download {0}',
    'viewer.copyseq': 'Copy the {0} sequence',
    'viewer.copied': 'copied',
    'viewer.copyfail': 'could not copy — select it in the table',
    'viewer.loading': 'Loading {0}…',
    'viewer.nostructure': 'This backbone has no structure file to show.',
    'viewer.noengine': 'Mol* did not load, so no structure can be drawn. Run '
      + 'python -m cyclicnano.webui.fetch_vendor on the server to fetch it locally.',
    'viewer.initfail': 'Mol* failed to start: {0}',
    'viewer.loadfail': 'Could not read the structure: {0}',
    'viewer.themefail': 'That colouring does not apply to this structure; '
      + 'Mol*’s default is shown instead.',
    'viewer.diagnose': 'Viewer self-check',
    'viewer.diagnose.hint': 'The result was also sent to the server log, so it can be '
      + 'read without copying anything out of the browser.',

    /* ----------------------------------------------------------- the jobs */
    'jobs.title': 'Submit a run',
    'jobs.sub': 'Fill the parameters in, click once, and the run starts on this '
      + 'machine; the log follows on the next page. The ranges are enforced by the '
      + 'server — the limits here are a convenience.',
    'jobs.disabled': 'Submission is not enabled on this service.',
    'jobs.disabled.why': 'It is off by default because the page has no '
      + 'authentication: anyone who can reach this port could spend GPU hours on a '
      + 'machine other people are using. Start it with --allow-submit and point '
      + '--profile at the machine configuration.',
    'jobs.params': 'Parameters',
    'jobs.preset': 'Filter preset',
    'jobs.preset.help': 'The stage-01 rule set, applied as one more config file.',
    'jobs.preset.none': 'no filter (every backbone continues)',
    'jobs.stages': 'Run through',
    'jobs.stages.help': 'Stages 00–01 alone take minutes and tell you how many '
      + 'backbones survive the filter. Finding that out after the whole pipeline has '
      + 'run costs hours.',
    'jobs.stages.all': 'everything, 00–05 (with interface scoring)',
    'jobs.stages.through_delivery': '00–04 (no interface scoring)',
    'jobs.stages.backbones_only': '00–01 only (see what the filter leaves)',
    'jobs.start': 'Start the run',
    'jobs.busy': 'A job is running, so a submit will be refused.',
    'jobs.free': 'One job at a time.',
    'jobs.history': 'Recent jobs',
    'jobs.warn.above': '{0} is set to {1}, above the usual {2}. This is the longest '
      + 'GPU stage, so be sure before starting.',
    'jobs.warn.full': 'The full pipeline includes ProteinMPNN, AlphaFold2 and a '
      + 'Rosetta relax. Measured at 500 backbones, that is a matter of hours.',
    'job.state.running': 'running',
    'job.state.finished': 'finished',
    'job.state.failed': 'failed',
    'job.state.cancelled': 'cancelled',
    'job.state.unknown': 'state unknown',
    'job.resumed': 'resuming an existing run',
    'job.id': 'Job',
    'job.started': 'Started',
    'job.elapsed': 'Elapsed',
    'job.stages': 'Stages',
    'job.returncode': 'Exit code',
    'job.log': 'Log',
    'job.cancel': 'Cancel the job',
    'job.follow': 'New output is pulled every 3 seconds.',
    'job.showlog': 'See the log →',
    'job.command': 'The full command (paste it into a terminal to repeat the run)',
    'job.title': 'Job {0}',
    'time.seconds': '{0} s',
    'time.minutes': '{0} min',
    'time.hours': '{0} h',
    'none': 'none',
  },

  zh: {
    'projects.title': '项目',
    'projects.sub': '一个 project 就是一个运行目录：设计目标、在它上面跑过的阶段、'
      + '以及这些阶段产出的全部结果。点一个看它交付的骨架。',
    'projects.new': '新建项目',
    'projects.none': '还没有项目',
    'projects.none.detail': '输出目录下没有像项目的子目录。一个目录只要有 '
      + 'project.json、tables/backbones.csv 或 04_delivery/delivery.json 之一就算。'
      + '下划线开头的目录不会显示。',
    'projects.empty': '还没有跑过任何阶段',
    'projects.empty.hint': '这个项目只有名字和目标，还没有结果。提交一个任务就能开始跑。',
    'projects.mismatch': '声明为 {0}，但 tables 里测到的是 {1}',
    'projects.created': '创建于 {0}',

    'new.title': '新建项目',
    'new.sub': '项目就是输出根目录下的一个文件夹。新建只写入名字和目标，不做别的：'
      + '给目标起名字和为它花 GPU 时间是两件事，所以这里不会跑任何阶段。',
    'new.name': '目录名',
    'new.name.help': '会变成 <输出根目录>/<名字>/，也会出现在页面 URL 里。'
      + '只能用字母、数字、下划线、点和连字符。',
    'new.label': '标题',
    'new.label.help': '显示在卡片上。自由文本，留空就用目录名。',
    'new.symmetry': '对称性',
    'new.symmetry.help': 'Cn 的 n。会记进项目，提交任务时作为默认值。',
    'new.monomer': '单体长度 (aa)',
    'new.monomer.help': '每条链的残基数。RFdiffusion 的 contig 总长会变成 n 倍于此。',
    'new.note': '备注',
    'new.note.help': '可选。这个项目是干什么的，用你自己的话写。',
    'new.create': '创建项目',
    'new.cancel': '返回',
    'new.disabled': '这个服务没有开启项目创建。',

    'jobs.field.name': '运行名称',
    'jobs.field.name.help': '结果写到 <输出根目录>/<名称>/。用已存在的名称会续跑那个运行，'
      + '已完成的阶段会跳过。',
    'jobs.field.symmetry': '对称性',
    'jobs.field.symmetry.help': 'Cn 的 n。contig 总长会变成 n × 单体长度。',
    'jobs.field.monomer_length': '单体长度 (aa)',
    'jobs.field.n_backbones': '生成骨架数',
    'jobs.field.n_backbones.help': '01 步的通过率 C5 约 2.5%、C3 约 9%，所以这个数要大。',
    'jobs.field.noise_scale': 'RFdiffusion noise',
    'jobs.field.noise_scale.help': '0.5 比 1.0 多出一倍的三段式螺旋；'
      + '0.0 螺旋更多但多是长发夹，已实测排除。',
    'jobs.field.batch_size': '每批骨架数',
    'jobs.field.batch_size.help': '分批只是为了断点续跑，批边界本身不可复现。',
    'jobs.field.n_seq_per_backbone': '每骨架序列数',
    'jobs.field.rmsd_cutoff': 'RMSD 阈值 (Å)',
    'jobs.field.rmsd_cutoff.help': '判定序列是否折回设计出的单体。文献常用 2.0，'
      + '这条流水线按要求用 1.0。',
    'jobs.field.plddt_cutoff': 'pLDDT 阈值',
    'jobs.field.seed': '随机种子',
    'jobs.submit.here': '提交任务',
    'jobs.submit.locked': '运行名称和对称性由这个项目决定。',

    'iface.title': '界面残基',
    'iface.sub': '共 {0} 个残基，按两条链间重原子 {1} Å 以内计',
    'iface.chain': '链 {0}',
    'iface.none': '两条链之间没有找到接触。',
    'iface.loading': '正在计算界面…',
    'iface.fail': '读取界面失败：{0}',

    'nav.submit': '提交任务',
    'theme.toggle': '切换明暗主题',
    'lang.toggle': '切换语言',
    'lang.name': 'EN',
    'root.unreadable': '（不可读）',

    'runs.complete': '已完成',
    'runs.incomplete': '未标记完成',
    'runs.undelivered': '未到交付阶段',
    'runs.undelivered.title': '有 tables/，还没有 04_delivery/',
    'runs.rate': '整体通过率 {0}%',
    'runs.interface': '界面评价 {0} 条序列',
    'runs.updated': '更新于 {0}',

    'funnel.generated': '生成骨架',
    'funnel.generated.title': 'RFdiffusion 生成的对称骨架总数',
    'funnel.filtered': '过滤通过',
    'funnel.filtered.title': '通过 01 步骨架过滤的骨架数',
    'funnel.sequences': 'MPNN 序列',
    'funnel.sequences.title': 'ProteinMPNN 为通过的骨架设计的序列总数',
    'funnel.predictions': 'AF2 预测',
    'funnel.predictions.title': '完成无 MSA 单体结构预测的序列数',
    'funnel.folded': '成功序列',
    'funnel.folded.title': 'RMSD 与 pLDDT 都达标的序列数',
    'funnel.delivered': '交付骨架',
    'funnel.delivered.title': '成功率达标、最终交付的骨架数',
    'funnel.kept': '，上一步的 {0}%',

    'run.sub': '{0} 个通过过滤的骨架，其中 {1} 个成功率达标{2}{3}。'
      + '点一行看它的序列和结构；表头可以点，用来排序。',
    'run.sub.monomer': ' · 单体 {0} aa',
    'run.nodelivery.title': '{0} · 还没有交付表',
    'run.nodelivery.sub': '这个运行还没跑到 04_delivery，所以没有可以逐骨架看的表。'
      + '每一步筛掉多少，在运行列表的漏斗里能看到。',
    'run.nodelivery.hint': '要补齐交付表，用同样的运行名称再提交一次、阶段选到 04 或 05 '
      + '即可；已经跑完的阶段会跳过，不会重复花 GPU。',

    'col.rank': '排名',
    'col.rank.title': '交付排名；成功率未达标的骨架没有排名',
    'col.backbone': '骨架',
    'col.delivered': '交付',
    'col.delivered.title': '单体成功率是否达标',
    'col.success': '成功率',
    'col.success.title': '该骨架的序列里 RMSD 与 pLDDT 都达标的比例',
    'col.rmsd.dist.title': '该骨架全部序列按与设计骨架的 CA RMSD 画成密度曲线，'
      + '横轴是所有行共用的 0–4 Å，曲线下方的细刻度是 20 个真实值的位置。'
      + '两端的数字是这个骨架自己的极小值和极大值，竖线是阈值。'
      + '高度按行各自归一化，所以峰高表示「在这里更集中」，不表示「序列更多」。',
    'col.plddt.dist.title': '该骨架全部序列按平均 pLDDT 画成密度曲线，'
      + '横轴是共用的 50–100，下方细刻度是真实值位置。'
      + '两端数字是这个骨架自己的极小值和极大值，竖线是阈值。',
    'col.ddg.dist.title': '已评价的序列按界面能画成密度曲线，横轴是共用的 −85–+5 REU，'
      + '下方细刻度是真实值位置。两端数字是这个骨架自己的极小值和极大值，竖线是阈值。'
      + '没评价过的序列不计入，数量在悬停提示里。',
    'col.iface.rate': '界面成功率',
    'col.iface.rate.title': '已评价的序列里 ddG 低于阈值的比例。'
      + '与单体成功率并列展示，不会混进去。',
    'col.seqs': '序列',
    'col.seqs.title': 'ProteinMPNN 设计的序列数',
    'col.preds': '预测',
    'col.preds.title': '完成 AF2 单体预测的序列数',
    'col.failreason': '未通过原因',
    'col.index': '#',
    'col.pass': '成功',
    'col.pass.title': 'RMSD 与 pLDDT 是否都达标',
    'col.rmsd': 'RMSD',
    'col.rmsd.title': '预测单体与设计骨架叠合后的 CA RMSD（Å）',
    'col.plddt': 'pLDDT',
    'col.plddt.title': '预测单体的平均 pLDDT',
    'col.ddg': 'ddG',
    'col.ddg.title': '列名按组里的习惯写作 ddG。这一列装的是 Rosetta '
      + 'InterfaceAnalyzerMover 的 dG_separated：复合物的能量减去两条链分开并重新 '
      + 'pack 之后的能量之和，单位 REU（不是 kcal/mol）。严格讲它是结合能 ΔG 而不是 '
      + 'ΔΔG——ΔΔG 是两个 ΔG 之差，需要一个参考态，de novo 设计没有可减的参考。'
      + '它只决定 interface_pass，不参与其他任何筛选。',
    'col.ddg.open': '点一下切到这个能量所测的 relax 后二聚体，并列出界面残基',
    'col.iface': '界面',
    'col.iface.title': '这条序列的 ddG 是否低于阈值。只描述界面，不参与筛选。',
    'col.sequence': '序列',

    'export.scope': '导出范围',
    'export.fasta': '下载 FASTA',
    'export.csv': '下载 CSV（Excel 可直接打开）',
    'export.scope.delivered': '交付的骨架 · 成功序列',
    'export.scope.delivered.help': '成功率达标的骨架里，RMSD 与 pLDDT 都通过的序列。'
      + '要下单合成就用这一档。',
    'export.scope.passed': '全部成功序列',
    'export.scope.passed.help': '再加上那些骨架整体成功率没达标、但自己折回来了的序列'
      + '——运气好的单条设计藏在这里。',
    'export.scope.all': '全部序列（含失败）',
    'export.scope.all.help': '交付目录里的每一条，失败的也在，留档用。',
    'export.note': '导出文件里界面能那一列叫 dG_interface，和表格里显示的 {0} 是同一个数；'
      + '文件用字段原名，这样脱离这个页面之后还能追回它的来源。',

    'detail.sequences': '序列',
    'detail.structure': '结构',
    'detail.criteria': '序列成功判定：RMSD ≤ {0} Å 且 平均 pLDDT ≥ {1}；'
      + '这两条决定成功率和交付判定。界面另算，{2} < {3} REU，不参与任何筛选。',
    'detail.criteria.noddg': '序列成功判定：RMSD ≤ {0} Å 且 平均 pLDDT ≥ {1}。'
      + '{2} 只描述界面，不参与筛选。',
    'fact.success': '成功率',
    'fact.passed': '成功序列',
    'fact.bestrmsd': '最佳 RMSD',
    'fact.bestplddt': '最佳 pLDDT',
    'fact.iface': '界面成功率',
    'fact.ifacecount': '已评价 {1} 条里 {0} 条',
    'fact.delivered': '交付',
    'yes': '是',
    'no': '否',
    'detail.pickrow': '点一行切换右边显示的结构。',
    'detail.best': '最佳',

    'viewer.structure': '结构',
    'viewer.colour': '着色',
    'viewer.group.backbone': '设计骨架',
    'viewer.group.current': '当前序列 {0}',
    'viewer.group.other': '其他序列',
    'viewer.missing': '（文件缺失）',
    'struct.backbone': '设计骨架',
    'struct.monomer': 'AF2 单体预测',
    'struct.oligomer': '{0} 聚体（设计几何）',
    'struct.oligomer.3': '三聚体（设计几何）',
    'struct.oligomer.5': '五聚体（设计几何）',
    'struct.dimer': 'relax 后的二聚体',
    'note.backbone': 'RFdiffusion 生成并通过过滤的骨架，只有主链，没有侧链。',
    'note.monomer': '无 MSA 的单体预测，B 因子列装的是 pLDDT。AlphaFold2 用自己的坐标系'
      + '输出，所以它不落在设计骨架上，RMSD 是叠合之后才算的。',
    'note.oligomer': '把单体预测按设计骨架的每个亚基位置摆一份，在这里拼成整体。'
      + '只做刚体叠合、没有 relax，相邻亚基之间可能有碰撞。流水线里没有任何一步'
      + '预测聚体：这是设计要求的几何，带上该序列对应的侧链，不构成该序列真能组装'
      + '成这个样子的证据。',
    'note.dimer': '把单体预测叠到设计骨架的 A、B 两条链上，再做带约束的 relax。'
      + '{0} 就是在这个结构上算的，所以它量的是设计出来的那个界面。relax 会把整体'
      + '平移旋转 2～3 Å，但内部几何只变 0.8 Å、两链间距变化不到 1 Å'
      + '——漂移是刚体的，界面几何守住了。下面列出的那些残基在 viewer 里是选中态，'
      + '点视口里的空白处可以取消。',
    'colour.chain': '按链',
    'colour.plddt': '按 pLDDT',
    'colour.ss': '按二级结构',
    'colour.na': '（此结构无此信息）',
    'legend.plddt': 'pLDDT',
    'legend.plddt.vhigh': '> 90 很高',
    'legend.plddt.high': '70–90 可信',
    'legend.plddt.low': '50–70 偏低',
    'legend.plddt.vlow': '≤ 50 很低',
    'legend.plddt.source': '（AlphaFold 官方配色，取自 B 因子列）',
    'legend.molstar': '配色由 Mol* 的 {0} 主题决定。',
    'viewer.download': '下载 {0}',
    'viewer.copyseq': '复制 {0} 的序列',
    'viewer.copied': '已复制',
    'viewer.copyfail': '复制失败，请从表格里选',
    'viewer.loading': '正在载入 {0} …',
    'viewer.nostructure': '这个骨架没有可显示的结构文件。',
    'viewer.noengine': 'Mol* 没有载入，结构无法显示。在服务器上运行 '
      + 'python -m cyclicnano.webui.fetch_vendor 把它取到本地。',
    'viewer.initfail': 'Mol* 初始化失败：{0}',
    'viewer.loadfail': '读取结构失败：{0}',
    'viewer.themefail': '这个结构用不了所选的着色，已退回 Mol* 的默认配色。',
    'viewer.diagnose': 'viewer 自检',
    'viewer.diagnose.hint': '自检结果同时发到了服务器日志，不用从浏览器里抄出来。',

    'jobs.title': '提交任务',
    'jobs.sub': '参数填好点一次，任务就在这台机器上跑起来，日志在下一页实时看。'
      + '范围由服务端校验，这里的限制只是方便。',
    'jobs.disabled': '这个服务没有开启任务提交。',
    'jobs.disabled.why': '它默认关闭，因为页面没有认证：能访问这个端口的人都能在这台'
      + '共用机器上起 GPU 任务。要开启就在启动命令里加 --allow-submit，'
      + '并用 --profile 指定机器配置。',
    'jobs.params': '参数',
    'jobs.preset': '骨架过滤预设',
    'jobs.preset.help': '01 步的筛选规则，作为一个额外的配置文件叠上去。',
    'jobs.preset.none': '不过滤（全部骨架进入下一步）',
    'jobs.stages': '跑到哪一步',
    'jobs.stages.help': '先只跑 00–01 只花几分钟就能看出过滤剩多少；'
      + '跑完整套之后才发现太少，代价是几个小时。',
    'jobs.stages.all': '全流程 00–05（含界面评价）',
    'jobs.stages.through_delivery': '00–04（不做界面评价）',
    'jobs.stages.backbones_only': '只到 00–01（先看过滤剩多少）',
    'jobs.start': '开始运行',
    'jobs.busy': '当前有任务在跑，提交会被拒绝。',
    'jobs.free': '一次只跑一个任务。',
    'jobs.history': '历史任务',
    'jobs.warn.above': '{0} 设成了 {1}，超过常用的 {2}。'
      + '这一步是 GPU 上最久的一段，先想清楚再跑。',
    'jobs.warn.full': '全流程含 ProteinMPNN、AlphaFold2 和 Rosetta relax；'
      + '按之前 500 骨架的实测，这是几个小时的量级。',
    'job.state.running': '运行中',
    'job.state.finished': '已完成',
    'job.state.failed': '失败',
    'job.state.cancelled': '已取消',
    'job.state.unknown': '状态未知',
    'job.resumed': '续跑已存在的运行',
    'job.id': '任务编号',
    'job.started': '开始',
    'job.elapsed': '用时',
    'job.stages': '阶段',
    'job.returncode': '退出码',
    'job.log': '日志',
    'job.cancel': '取消任务',
    'job.follow': '日志每 3 秒拉取一次新增的部分。',
    'job.showlog': '看日志 →',
    'job.command': '完整命令（可以直接贴到终端重跑）',
    'job.title': '任务 {0}',
    'time.seconds': '{0} 秒',
    'time.minutes': '{0} 分钟',
    'time.hours': '{0} 小时',
    'none': '无',
  },
};

const LANG_KEY = 'cyclicnano.lang';

/** The chosen language, English unless the viewer asked for the other one. */
let lang = 'en';
try {
  const stored = localStorage.getItem(LANG_KEY);
  if (stored && STRINGS[stored]) lang = stored;
} catch (e) { /* private window, or site data blocked */ }

function setLang(next) {
  if (!STRINGS[next]) return;
  lang = next;
  document.documentElement.lang = next === 'zh' ? 'zh-CN' : 'en';
  try { localStorage.setItem(LANG_KEY, next); } catch (e) { /* blocked */ }
}

const otherLang = () => (lang === 'en' ? 'zh' : 'en');

/** Whether a string exists at all, in either table.
 *
 *  `t()` falls back to the key itself, which is the right answer for a label
 *  that someone forgot to add and the wrong one for an optional string: a job
 *  field with no help text would otherwise render `jobs.field.seed.help` at the
 *  reader. Callers of optional strings ask first.
 */
function has(key) {
  return key in (STRINGS[lang] || {}) || key in STRINGS.en;
}

/** One string, with `{0}`-style substitutions.
 *
 *  A key the chosen language is missing falls back to English, and a key neither
 *  has comes back as itself: a visible gap rather than an empty element.
 */
function t(key, ...args) {
  const table = STRINGS[lang] || STRINGS.en;
  const raw = (key in table) ? table[key] : (STRINGS.en[key] ?? key);
  return String(raw).replace(/\{(\d+)\}/g, (match, i) => {
    const value = args[Number(i)];
    return value === undefined || value === null ? match : String(value);
  });
}

document.documentElement.lang = lang === 'zh' ? 'zh-CN' : 'en';
