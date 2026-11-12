#!/usr/bin/env python3
"""snk-chipseq: run plotHeatmap or plotProfile with the profile legend placed
outside the plot (right of the axes). deepTools' --legendLocation only places it
inside the axes; this patches Axes.legend and pyplot.savefig, then calls the
tool's main() with the same CLI arguments.

Usage: deeptools_legend_outside.py {plotHeatmap|plotProfile} <tool arguments>
"""
import importlib
import sys

import matplotlib.axes
import matplotlib.pyplot as plt

_legend = matplotlib.axes.Axes.legend
_savefig = plt.savefig
_legends = []


def legend_outside(self, *args, **kwargs):
    kwargs.pop("loc", None)
    kwargs.pop("bbox_to_anchor", None)
    kwargs.update(loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0)
    leg = _legend(self, *args, **kwargs)
    # Keep constrained_layout / tight_layout from shrinking the panels for it.
    leg.set_in_layout(False)
    _legends.append(leg)
    return leg


def savefig(*args, **kwargs):
    # ...but grow the saved image to include it (plotProfile saves without
    # bbox_inches; plotHeatmap already uses bbox_inches='tight', pad 0.1).
    kwargs["bbox_inches"] = "tight"
    kwargs.setdefault("pad_inches", 0.1)
    kwargs.setdefault("bbox_extra_artists", _legends)
    return _savefig(*args, **kwargs)


matplotlib.axes.Axes.legend = legend_outside
plt.savefig = savefig

if __name__ == "__main__":
    tools = {"plotHeatmap": "deeptools.plotHeatmap", "plotProfile": "deeptools.plotProfile"}
    if len(sys.argv) < 2 or sys.argv[1] not in tools:
        sys.exit(f"usage: {sys.argv[0]} {{{'|'.join(tools)}}} <arguments>")
    importlib.import_module(tools[sys.argv[1]]).main(sys.argv[2:])
