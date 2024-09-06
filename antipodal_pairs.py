import argparse
from timeit import default_timer as timer
from datetime import timedelta
from string import ascii_uppercase

import matplotlib.pyplot as plt

from convexcpp import ROI


COORDS = {
    'box': [
        (1, 1),
        (5, 2),
        (5, 6),
        (1, 5),
    ],
    'pentagon': [
        (-3, -1),
        (3, -1),
        (3, 2),
        (0, 4),
        (-3, 2),
    ],
    'hexagon': [
        (-3, -1),
        (3, -1),
        (3, 2),
        (0, 4),
        (-3, 2),
        (0, -3),
    ],
    'shamos': [
        (-0.25, 0),
        (4, 2),
        (4, 7),
        (-1, 9),
        (-1.5, 6),
    ],
    'aoi': [
        (0.998055, 0.0),
        (1.258611, 0.226944),
        (0.260833, 1.213611),
        (0.0, 0.985),
    ]
}

# Parse command line arguments
parser = argparse.ArgumentParser()

parser.add_argument(
    'geom',
    type=str,
    help=f'ROI geometry. Valid: {", ".join(list(COORDS.keys()))}',
)

parser.add_argument(
    '-p',
    '--plot',
    dest='do_plot',
    action='store_true',
    help='Plot the ROI',
)

args = parser.parse_args()

# Compute runtime of ROI instantiation and antipodal pair computation
start_time = timer()

r = ROI(COORDS[args.geom])
pairs = r.get_antipodal_pairs()

end_time = timer()

micro_s = timedelta(seconds=end_time - start_time).microseconds
print(f'Execution time: {micro_s} {chr(956)}s ({micro_s * 0.001:.3f} ms)')

if True:
    for p in pairs:
        print(p)

# Get the coordinates of the paired vertices
antipodal_coords = [(r.vertices[i].coords, r.vertices[j].coords)
                    for i, j in pairs]

# Plot ROI outline, label ROI vertices, and connect antipodal pairs
if args.do_plot:
    fig = plt.figure()
    ax = fig.add_subplot(111)

    colors = 'bgrcmybgrcmybgrcmybgrcmy'

    v_coords = r.vertices_coords + [r.vertices_coords[0]]
    x, y = list(zip(*v_coords))
    ax.plot(x, y, 'black', zorder=0)

    for i, pt in enumerate(zip(x[:-1], y[:-1])):
        ax.text(*pt, ascii_uppercase[i], fontsize=14)

    for i, x in enumerate(antipodal_coords):
        ax.plot(*list(zip(*x)), color=colors[i])

    ax.grid()
    plt.show()
