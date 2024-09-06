from .util import is_close, Point


class Vertex:
    """A vertex of the Region of Interest."""

    def __init__(self, x: float, y: float, name: str) -> None:
        """Constructor.

        Parameters
        ----------
        x: float
            X-coordinate.
        y: float
            Y-coordinate.
        name: str
            Vertex name.
        """
        self.x = x
        self.y = y
        self.name = name

    @property
    def coords(self) -> Point:
        return self.x, self.y

    def __eq__(self, value: object) -> bool:
        return (
            isinstance(value, self.__class__) and
            value.name == self.name and
            is_close(value.x, self.x) and
            is_close(value.y, self.y)
        )

    def __repr__(self) -> str:
        """Return a string representation of the Vertex."""
        name = f'{self.name} ' if self.name != '' else self.name
        return f'<Vertex {name}({self.x:.3f}, {self.y:.3f})>'
