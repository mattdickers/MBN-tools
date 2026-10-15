import pandas as pd
import plotext as plt
import subprocess
import os
import argparse

# -------------------------
# Command-line arguments
# -------------------------

parser = argparse.ArgumentParser(
    description="Live plot data from an MBN Explorer chemical statistics file."
)

parser.add_argument(
    "filename",
    help="Path to the MBN Explorer .dat file"
)

parser.add_argument(
    "-x", "--x",
    default='Step',
    help="Column to use for the x-axis"
)

parser.add_argument(
    "-y", "--y",
    help="Column to use for the y-axis. If not specificed all columns will be plotted."
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

num_lines = subprocess.run(["wc", "-l", args.filename], capture_output=True, text=True, check=True)
num_lines = int(num_lines.stdout.split()[0])-1

n_tail = min(args.lines, num_lines)

fig = plt.figure
frame = 0
rows = []
checked = False

while True:
    w, h = plt.terminal.size(update=True)
    if frame:
        plt.terminal.clean(h)

    # Get updated tail of chem stats file
    file = subprocess.run(['tail', '-n', str(n_tail), args.filename], capture_output=True, text=True, check=True).stdout
    for line in file.split('\n')[:-1]:
        cols = line.split()

        data = {'Step': int(cols[0])}

        for i in range(1, len(cols), 2):
            particle_type = cols[i]
            count = int(cols[i+1])
            data[particle_type] = count
        rows.append(data)

    data = pd.DataFrame(rows)

    # Set up data for plotting
    x = data[args.x]

    # Plot data
    fig.clear.data()
    fig.plot_size(w, h-1)

    fig.theme(args.theme)

    # Axes labels
    fig.title('Press q to exit')
    fig.label(args.x, axis=0)
    fig.label(args.y, axis=1)

    fig.ruler("both", "both").grid(True)

    if args.y:
        # Check that requested columns exist
        if not checked:
            if args.x not in data.columns:
                parser.error(
                    f"X column '{args.x}' not found. "
                    f"Available columns: {', '.join(data.columns)}"
                )

            if args.y not in data.columns:
                parser.error(
                    f"Y column '{args.y}' not found. "
                    f"Available columns: {', '.join(data.columns)}"
                )

            checked = True
        
        y = data[args.y]
        fig.draw(fig.signal(x, y, marker='braille').lines())
    else:
        for column in data.columns:
            if column == 'Step':
                continue
            fig.draw(fig.signal(x, data[column], marker='braille').lines().label(column))

    fig.show(flush=True)

    if plt.terminal.is_pressed('q'):
        os.system('clear')
        break

    plt.sleep(args.interval)
    frame+=1
