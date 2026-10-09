import pandas as pd
import plotext as plt
import subprocess
from io import StringIO
import time
import os
import argparse

# -------------------------
# Command-line arguments
# -------------------------

parser = argparse.ArgumentParser(
    description="Live plot data from an MBN Explorer log file."
)

parser.add_argument(
    "filename",
    help="Path to the MBN Explorer .log file"
)

parser.add_argument(
    "--x",
    required=True,
    help="Column to use for the x-axis"
)

parser.add_argument(
    "--y",
    required=True,
    help="Column to use for the y-axis"
)

parser.add_argument(
    "-n", "--lines",
    type=int,
    default=50,
    help="Number of lines to plot (default: 10)"
)

parser.add_argument(
    "-i", "--interval",
    type=float,
    default=1.0,
    help="Update interval in seconds (default: 1)"
)

parser.add_argument(
    "-t", "--theme",
    type=str,
    default='default',
    help="plotext theme. For available options see: https://plotext.readthedocs.io/en/latest/theme.html"
)


args = parser.parse_args()


# Get log file headings
headings = subprocess.run(['head', '-n', '1', args.filename], text=True, capture_output=True).stdout.split()

# Check that requested columns exist
if args.x not in headings:
    parser.error(
        f"X column '{args.x}' not found. "
        f"Available columns: {', '.join(headings)}"
    )

if args.y not in headings:
    parser.error(
        f"Y column '{args.y}' not found. "
        f"Available columns: {', '.join(headings)}"
    )

num_lines = subprocess.run(["wc", "-l", args.filename], capture_output=True, text=True, check=True)
num_lines = int(num_lines.stdout.split()[0])-1

n_tail = min(args.lines, num_lines)

fig = plt.figure
frame = 0

while True:
    w, h = plt.terminal.size(update=True)
    if frame:
        plt.terminal.clean(h)

    # Get updated tail of log file
    log = subprocess.run(['tail', '-n', str(n_tail), args.filename], capture_output=True, text=True, check=True)
    data = pd.read_csv(StringIO(log.stdout), sep=r'\s+', header=None, names=headings)

    # Set up data for plotting
    x = data[args.x]
    y = data[args.y]

    # Plot data
    fig.clear.data()
    fig.plot_size(w, h - 1)

    fig.theme(args.theme) 

    # Axes labels
    fig.title('Press q to exit')
    fig.label(args.x, axis=0)
    fig.label(args.y, axis=1)

    fig.ruler("both", "both").grid(True)
    fig.draw(fig.signal(x, y).lines())
    fig.show(flush=True)

    if plt.terminal.is_pressed('q'):
        os.system('clear')
        break

    plt.sleep(args.interval)
    frame+=1
