"""Pipeline state shared by every tab.

The pipeline runs in one direction,

    DATA -> FILTER -> SEGMENT -> UNWRAP -> MOTION -> STRAIN

so invalidation is a sweep over the later stages: change a stage's settings
and every result after it is dropped. Tabs watch the signals and grey
themselves out, so a result can never be shown against stale inputs.
"""

from __future__ import annotations

from enum import IntEnum

from PySide6.QtCore import QObject, Signal

from .snapshots import SnapshotStore

__all__ = ["Stage", "Session"]


class Stage(IntEnum):
    DATA = 0
    FILTER = 1
    SEGMENT = 2
    UNWRAP = 3
    MOTION = 4
    STRAIN = 5

    @property
    def title(self) -> str:
        return {
            Stage.DATA: "Data",
            Stage.FILTER: "Filtering",
            Stage.SEGMENT: "Segmentation",
            Stage.UNWRAP: "Unwrapping",
            Stage.MOTION: "Motion",
            Stage.STRAIN: "Strain",
        }[self]


class Session(QObject):
    """Holds one configuration and one result per stage."""

    stageChanged = Signal(int)       # a stage's result appeared or vanished
    busyChanged = Signal(bool)
    frameRequested = Signal(int)     # a tab asks to show a particular frame
    message = Signal(str)
    snapshotsChanged = Signal()      # a snapshot was kept or dropped

    def __init__(self, parent=None):
        super().__init__(parent)
        self._results: dict[Stage, object] = {}
        self._configs: dict[Stage, object] = {}
        self._busy = False

        # Segmentation state, keyed by (slice, frame)
        self.contours: dict[tuple[int, int], list] = {}
        self.labels: dict[tuple[int, int], list] = {}
        self.mask = None

        # Where a pixel is sampled, in the canvas coordinates contours are drawn
        # in. The canvas puts pixel i over [i, i+1), so its centre is at i+0.5
        # and 0.5 is the value that makes a mask line up with what was drawn.
        # Projects written before this was fixed are pinned to 0.0 on load, so
        # reopening one reproduces the mask it was saved with.
        self.mask_origin = 0.5

        self.structure_masks: dict[str, object] = {}
        """Per-structure masks, keyed 'lv', 'rv' or 'other'. Analysis only."""
        self.mask_source = "none"
        """Where the current mask came from: 'none', 'contours' or 'file'."""
        self.dataset = None

        # Unwrapping settings, shared so the preview and the solver agree.
        # Previously the Unwrapping tab had controls that reached its own
        # display and nothing else.
        self.unwrap = {"connectivity": 4, "seed": "auto"}

        # Outside _results on purpose: invalidate_from must not reach these.
        # A snapshot exists to outlive the settings that produced it.
        self.snapshots = SnapshotStore()

        self.current_frame = 0
        self.current_slice = 0
        self.project_path = None
        self.dirty = False

    # ------------------------------------------------------------- results
    def result(self, stage: Stage):
        return self._results.get(stage)

    def has(self, stage: Stage) -> bool:
        return self._results.get(stage) is not None

    def set_result(self, stage: Stage, value):
        self._results[stage] = value
        self.invalidate_after(stage)
        self.dirty = True
        self.stageChanged.emit(int(stage))

    def config(self, stage: Stage):
        return self._configs.get(stage)

    def set_config(self, stage: Stage, cfg):
        if self._configs.get(stage) == cfg:
            return
        self._configs[stage] = cfg
        self.invalidate_from(stage)

    # -------------------------------------------------------- invalidation
    def invalidate_from(self, stage: Stage):
        """Drop this stage's result and every later one."""
        changed = False
        for s in Stage:
            if s >= stage and self._results.get(s) is not None:
                self._results[s] = None
                changed = True
        if changed:
            self.dirty = True
            self.stageChanged.emit(int(stage))

    def invalidate_after(self, stage: Stage):
        """Drop every result after this stage."""
        changed = False
        for s in Stage:
            if s > stage and self._results.get(s) is not None:
                self._results[s] = None
                changed = True
        if changed:
            self.stageChanged.emit(int(stage) + 1)

    def ready_for(self, stage: Stage) -> bool:
        """Whether every stage before this one has a result."""
        if stage == Stage.DATA:
            return True
        if stage == Stage.SEGMENT:
            return self.has(Stage.DATA)
        # Unwrapping is advisory. The solver unwraps internally, so requiring
        # the preview first only forced the same work to be done twice.
        advisory = {Stage.SEGMENT, Stage.UNWRAP}
        return all(self.has(s) for s in Stage if s < stage and s not in advisory) and (
            self.mask is not None or stage <= Stage.SEGMENT
        )

    # --------------------------------------------------------------- busy
    @property
    def busy(self) -> bool:
        return self._busy

    def set_busy(self, value: bool):
        if self._busy != value:
            self._busy = bool(value)
            self.busyChanged.emit(self._busy)
