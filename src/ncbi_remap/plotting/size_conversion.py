from typing import Tuple, Union

from pint import UnitRegistry


class FigSize:
    """Does simple figure size unit conversions."""
    ureg = UnitRegistry()

    def __init__(self, figsize: Tuple[float, float]):
        self._width = float(figsize[0]) * self.ureg.inch
        self._height = float(figsize[1]) * self.ureg.inch
        self._width_cm = self._width.to(self.ureg.cm)
        self._height_cm = self._height.to(self.ureg.cm)

    @property
    def width_in(self, as_string=True) -> Union[str, float]:
        if as_string:
            return f"{self._width.magnitude}in"
        return self._width.magnitude

    @property
    def height_in(self, as_string=True) -> Union[str, float]:
        if as_string:
            return f"{self._height.magnitude}in"
        return self._height.magnitude

    @property
    def width_cm(self, as_string=True) -> Union[str, float]:
        if as_string:
            return f"{self._width_cm.magnitude}cm"
        return self._width_cm.magnitude

    @property
    def height_cm(self, as_string=True) -> Union[str, float]:
        if as_string:
            return f"{self._height_cm.magnitude}cm"
        return self._height_cm.magnitude

