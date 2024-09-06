from __future__ import annotations

from typing import Tuple, TYPE_CHECKING

import numpy as np
from shapely.geometry import LineString

from .util import distance, is_close, Point

if TYPE_CHECKING:
    from shapely.geometry.base import BaseGeometry
    from .vertex import Vertex


class Edge:
    """An edge of the Region of Interest."""

    def __init__(self, origin: Vertex, terminus: Vertex, name: str) -> None:
        """Constructor.

        Parameters
        ----------
        origin: Tuple[float, float]
            First Edge vertex.
        terminus: Tuple[float, float]
            Second Edge vertex.
        name: str
            Edge name.
        """
        self.origin = origin
        self.terminus = terminus
        self.name = name
        self.linestring = LineString([origin.coords, terminus.coords])
        self.length = distance(origin.coords, terminus.coords)
        self.azimuth = self._calc_azimuth()

    @property
    def back_azimuth(self) -> float:
        """Return the back azimuth (azimuth + 180)."""
        return (self.azimuth + 180) % 360.

    @property
    def coords(self) -> Tuple[Point, Point]:
        """Return the coordinates of the origin and terminus vertices."""
        return self.origin.coords, self.terminus.coords

    @property
    def vertices(self) -> Tuple[Point]:
        """Return the vertices as a tuple."""
        return self.origin, self.terminus

    def _calc_azimuth(self) -> float:
        """Calculate the azimuth from point 1 to point 2."""
        x1, y1 = self.origin.coords
        x2, y2 = self.terminus.coords

        theta = np.arctan2(x2 - x1, y2 - y1)

        if theta < 0:
            theta += (2 * np.pi)

        return round(np.degrees(theta))

    def angle_between(
        self,
        edge: Edge,
        how: str = 'counter-clockwise',
    ) -> float:
        """Compute the angle between the azimuth of this Edge and
        another Edge.

        Parameters
        ----------
        edge: Edge
            The other edge.
        how: str, optional
            The direction to compute the azimuth difference.
            Either 'counter-clockwise' (default) or 'clockwise'.

        Returns
        -------
        float
        """
        if how == 'counter-clockwise':
            return (self.azimuth - edge.azimuth) % 360.
        elif how == 'clockwise':
            return (edge.azimuth - self.azimuth) % 360.
        else:
            raise ValueError(f'Invalid `how` argument "{how}"')

    def crosses(self, geom: BaseGeometry) -> bool:
        """Whether the Edge crosses another geometric object.

        From Shapely docs:
        "...if the interior of the object intersects the interior of the other
        but does not contain it, and the dimension of the intersection is less
        than the dimension of the one or the other"

        https://shapely.readthedocs.io/en/stable/manual.html#object.crosses

        Parameters
        ----------
        geom: shapely.Geometry
            A geometric object.

        Returns
        -------
        bool
            Whether the Edge crosses another geometric object.
        """
        return self.linestring.crosses(geom)

    def intersects(self, geom: BaseGeometry) -> bool:
        """Whether the Edge intersects another geometric object.

        From Shapely docs:
        "Geometric objects intersect if they have any boundary or
        interior point in common."

        https://shapely.readthedocs.io/en/stable/manual.html#object.intersects

        Parameters
        ----------
        geom: shapely.Geometry
            A geometric object.

        Returns
        -------
        bool
            Whether the Edge intersects another geometric object.
        """
        return self.linestring.intersects(geom)

    def intersection(self, geom: BaseGeometry) -> BaseGeometry:
        """Return the intersection of the Edge with another geometric object.

        Parameters
        ----------
        geom: shapely.Geometry
            A geometric object.

        Returns
        -------
        shapely.BaseGeometry
            The intersection of the edge's linestring and the given geometry.
        """
        return self.linestring.intersection(geom)

    def __eq__(self, value: object) -> bool:
        return (
            isinstance(value, self.__class__) and
            self.name == value.name and
            is_close(self.point_1[0], value.point_1[0]) and
            is_close(self.point_1[1], value.point_1[1]) and
            is_close(self.point_2[0], value.point_2[0]) and
            is_close(self.point_2[1], value.point_2[1])
        )

    def __repr__(self) -> str:
        """Return a string representation of the Edge."""
        def round_vtx(vtx):
            return round(vtx[0]), round(vtx[1])

        coords = [round_vtx(v.coords) for v in self.vertices]

        return f'<Edge "{self.name}" ({", ".join(map(str, coords))})>'
