from __future__ import annotations
from string import ascii_uppercase
from typing import Dict, List, Tuple, TYPE_CHECKING

import numpy as np
from scipy.spatial import ConvexHull
from shapely.geometry import Polygon

from .edge import Edge
from .vertex import Vertex
from .util import next_index, Point

if TYPE_CHECKING:
    from shapely.geometry.base import BaseGeometry


class ROI:
    """Convex Region of Interest."""

    def __init__(self, vertices: List[Point]) -> None:
        """ROI Constructor.

        Parameters
        ----------
        vertices: List of Tuple[float, float]
            List of tuples defining the x- and y-coordinates of the ROI
            vertices.
        """
        # Use ConvexHull to ensure vertices are in counter-clockwise order
        hull = ConvexHull(vertices)
        self._polygon = Polygon(hull.points[hull.vertices])
        self._vertices = self._init_vertices()
        self._edges = self._init_edges()

    @property
    def polygon(self) -> Polygon:
        """Return the Polygon."""
        return self._polygon

    @property
    def vertices(self) -> Dict[str, Vertex]:
        """Return the dict of Vertex objects."""
        return self._vertices

    @property
    def vertices_coords(self) -> List[Point]:
        """Return a list of Vertex coordinate tuples."""
        return [x.coords for x in self._vertices.values()]

    @property
    def edges(self) -> Dict[str, Edge]:
        """Return the dict of Edge objects."""
        return self._edges

    @property
    def num_edges(self) -> int:
        """Return the number of edges."""
        return len(self._edges.keys())

    @property
    def num_vertices(self) -> int:
        """Return the number of vertices."""
        return len(self._vertices.keys())

    def _init_vertices(self) -> List[Vertex]:
        """Initialize the ROI Vertices."""
        vertices = {}
        for idx, v in enumerate(list(self._polygon.exterior.coords)[:-1]):
            v_name = ascii_uppercase[idx]
            vertices[v_name] = Vertex(*v, v_name)

        return vertices

    def _init_edges(self) -> List[Edge]:
        """Initialize the ROI Edges."""
        edges = {}

        vertices = self.vertices.values()
        max_idx = len(vertices)

        for n in range(max_idx):
            vert1 = self.vertices[ascii_uppercase[n]]
            vert2 = self.vertices[ascii_uppercase[(n + 1) % max_idx]]
            e_name = f'{vert1.name}{vert2.name}'

            edges[e_name] = Edge(vert1, vert2, e_name)

        return edges

    def get_antipodal_pairs(self) -> List[Tuple[str, str]]:
        """Compute the antipodal pairs for all vertices.

        Parameters
        ----------
        None.

        Returns
        -------
        list of tuple of str, str
            Antipodal pairs identified by vertex names.
        """
        pairs = []
        n = self.num_vertices

        if n < 3:
            return []

        vertices = list(self.vertices.keys())

        i = 0
        j = 1

        edge_i = self.get_edge_from_vtx(vertices[i])
        edge_j = self.get_edge_from_vtx(vertices[j])

        while edge_i.angle_between(edge_j) < 180.:
            j = next_index(j, n)
            edge_j = self.get_edge_from_vtx(vertices[j])

        pairs.append((vertices[i], vertices[j]))

        while j != 1:
            a = 360. - (self.get_edge_from_vtx(vertices[i])
                        .angle_between(self.get_edge_from_vtx(vertices[j])))

            if a == 180:
                pairs.extend([
                    (vertices[next_index(i, n)], vertices[j]),
                    (vertices[i], vertices[next_index(j, n)]),
                    (vertices[next_index(i, n)], vertices[next_index(j, n)]),
                ])
                i = next_index(i, n)
                j = next_index(j, n)
            elif a < 180:
                pairs.append((vertices[next_index(i, n)], vertices[j]))
                i = next_index(i, n)
            else:
                pairs.append((vertices[i], vertices[next_index(j, n)]))
                j = next_index(j, n)

        return pairs

    def get_edge_from_vtx(self, vertex: str, point: str = 'origin') -> Edge:
        """Return an Edge given the name of its origin or terminus vertex.

        Parameters
        ----------
        vertex: str
            Name of the origin or terminus vertex.
        point: str, optional
            Which edge vertex to use for name matching.
            Either 'origin' or 'terminus'. Default is 'origin'.

        Returns
        -------
        Edge
        """
        if point == 'origin':
            return [x for x in self.edges.values() if x.name[0] == vertex][0]
        elif point == 'terminus':
            return [x for x in self.edges.values() if x.name[1] == vertex][0]
        else:
            raise ValueError(f'Invalid point argument "{point}"')

    def crosses(self, geom: BaseGeometry) -> bool:
        """Whether the ROI crosses another geometric object.

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
            Whether the ROI crosses another geometric object.
        """
        return self.linestring.crosses(geom)

    def intersects(self, geom: BaseGeometry) -> bool:
        """Whether the ROI intersects another geometric object.

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
            Whether the ROI intersects another geometric object.
        """
        return self.linestring.intersects(geom)

    def intersection(self, geom: BaseGeometry) -> BaseGeometry:
        """Return the intersection of the ROI with another geometric object.

        Parameters
        ----------
        geom: shapely.Geometry
            A geometric object.

        Returns
        -------
        shapely.BaseGeometry
            The intersection of the ROI's linestring and the given geometry.
        """
        return self.linestring.intersection(geom)

    def __repr__(self) -> str:
        """Return a string representation of the Grid instance."""
        coords = [(round(x, 3), round(y, 3))
                  for x, y in list(self._polygon.exterior.coords)][:-1]
        return f'<ROI ({", ".join(map(str, coords))})>'
