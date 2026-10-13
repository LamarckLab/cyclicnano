"""Declarative filter engine.

Filters live in YAML as plain expressions over the metric names, e.g.

    backbone:
      - "n_clash == 0"
      - "helix_frac >= 0.40"
      - "loop_max_len <= 12"

Nothing is ever deleted: every candidate keeps every metric and gains a pass flag plus
the list of rules it failed. Re-tuning a threshold therefore only replays this stage,
never the GPU stages above it.
"""
from __future__ import annotations

import ast
import operator
from dataclasses import dataclass

_BINOPS = {
    ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul,
    ast.Div: operator.truediv, ast.Pow: operator.pow, ast.Mod: operator.mod,
}
_CMPOPS = {
    ast.Lt: operator.lt, ast.LtE: operator.le, ast.Gt: operator.gt,
    ast.GtE: operator.ge, ast.Eq: operator.eq, ast.NotEq: operator.ne,
}
_FUNCS = {"abs": abs, "min": min, "max": max, "len": len, "round": round}


class RuleError(ValueError):
    """Raised when a rule is malformed or references an unknown metric."""


def _eval(node: ast.AST, ctx: dict):
    if isinstance(node, ast.Expression):
        return _eval(node.body, ctx)
    if isinstance(node, ast.Constant):
        return node.value
    if isinstance(node, ast.Name):
        if node.id not in ctx:
            raise RuleError(f"unknown metric: {node.id}")
        return ctx[node.id]
    if isinstance(node, ast.UnaryOp):
        val = _eval(node.operand, ctx)
        if isinstance(node.op, ast.Not):
            return not val
        if isinstance(node.op, ast.USub):
            return -val
        if isinstance(node.op, ast.UAdd):
            return +val
        raise RuleError(f"unsupported unary op: {type(node.op).__name__}")
    if isinstance(node, ast.BinOp):
        fn = _BINOPS.get(type(node.op))
        if fn is None:
            raise RuleError(f"unsupported operator: {type(node.op).__name__}")
        return fn(_eval(node.left, ctx), _eval(node.right, ctx))
    if isinstance(node, ast.BoolOp):
        vals = (_eval(v, ctx) for v in node.values)
        return all(vals) if isinstance(node.op, ast.And) else any(vals)
    if isinstance(node, ast.Compare):
        left = _eval(node.left, ctx)
        for op, right_node in zip(node.ops, node.comparators):
            fn = _CMPOPS.get(type(op))
            if fn is None:
                raise RuleError(f"unsupported comparison: {type(op).__name__}")
            right = _eval(right_node, ctx)
            if not fn(left, right):
                return False
            left = right
        return True
    if isinstance(node, ast.Call):
        if not isinstance(node.func, ast.Name) or node.func.id not in _FUNCS:
            raise RuleError("only abs/min/max/len/round may be called in a rule")
        return _FUNCS[node.func.id](*[_eval(a, ctx) for a in node.args])
    raise RuleError(f"unsupported syntax in rule: {type(node).__name__}")


def evaluate(rule: str, metrics: dict) -> bool:
    """Evaluate one rule against a metric dict. No builtins, no attribute access."""
    try:
        tree = ast.parse(rule, mode="eval")
    except SyntaxError as exc:
        raise RuleError(f"cannot parse rule {rule!r}: {exc}") from exc
    result = _eval(tree, metrics)
    if not isinstance(result, bool):
        raise RuleError(f"rule {rule!r} returned {type(result).__name__}, expected a boolean")
    return result


@dataclass
class FilterResult:
    passed: bool
    failed_rules: list[str]

    @property
    def reason(self) -> str:
        return "; ".join(self.failed_rules)


class FilterSet:
    """An ordered, named set of rules applied to one candidate at a time."""

    def __init__(self, name: str, rules: list[str]):
        self.name = name
        self.rules = list(rules or [])

    def check_syntax(self, sample_metrics: dict) -> None:
        """Fail fast at startup rather than after hours of GPU time."""
        for rule in self.rules:
            evaluate(rule, sample_metrics)

    def apply(self, metrics: dict) -> FilterResult:
        failed = [r for r in self.rules if not evaluate(r, metrics)]
        return FilterResult(passed=not failed, failed_rules=failed)

    def summary(self, rows: list[dict]) -> dict:
        """Per-rule rejection counts, so it is visible which threshold is doing the work."""
        counts = {r: 0 for r in self.rules}
        for m in rows:
            for rule in self.rules:
                if not evaluate(rule, m):
                    counts[rule] += 1
        n = len(rows)
        return {
            # The rules are repeated verbatim because the summary is the artefact that
            # ends up beside the results: a reader should not have to find the config
            # that produced a run to know what it selected on.
            "rules": list(self.rules),
            "n_input": n,
            "n_pass": sum(1 for m in rows if self.apply(m).passed),
            "rejected_by_rule": counts,
        }
