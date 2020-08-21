import numpy as np
import matplotlib.pyplot as plt
import copy
from matplotlib.patches import Rectangle
import os


def check_length(s):
    if len(s) == 81:
        return True
    else:
        return False


def unlist(l):
    return [item for sublist in l for item in sublist]


def str_to_array(string):
    s_rows = [string[i:i+9] for i in range(0, 81, 9)]
    s_split = [[int(y) for y in x] for x in s_rows]
    return np.array(s_split, dtype=np.uint8)


def get_grid_matrix(grid_str):
    if not check_length(grid_str):
        raise ValueError('Invalid grid length')
    for n in range(9):
        if grid_str.count(str(n)) != 9:
            raise ValueError('Invalid grid')
    return str_to_array(grid_str)


class Sudoku:
    def __init__(self, sudoku_string, grid_string=None, diag=False):
        if grid_string is None:
            grid_string = '000111222000111222000111222333444555333444555333444555666777888666777888666777888'
        self.grid = get_grid_matrix(grid_string)
        self.sudoku = str_to_array(sudoku_string)
        self.diag = diag
        self.guess_coords = []

    def check_cell(self, x, y, options=None, other=None, other_xy=None):
        if self.sudoku[x, y] == 0:
            if options is None:
                options = list(range(1, 10))
            if other_xy is not None:
                options = other.check_cell(other_xy[0], other_xy[1])
            options = [i for i in options if i not in list(self.sudoku[x, :])]
            options = [i for i in options if i not in list(self.sudoku[:, y])]
            options = [i for i in options if i not in list(self.sudoku[self.grid == self.grid[x, y]].flatten())]
            if self.diag and (x == y):
                options = [i for i in options if i not in [self.sudoku[j, j] for j in range(9)]]
            elif self.diag and ((x + y) == 8):
                options = [i for i in options if i not in [self.sudoku[j, abs(8-j)] for j in range(9)]]
            return options
        else:
            return []

    def check_grid(self, grid_nr, other=None, overlapping_coords=None, plot_progress=False):
        if not np.any(self.sudoku[self.grid == grid_nr] == 0):
            grid_coords = [tuple(x) for x in list(np.array(np.where(self.grid == grid_nr)).transpose())]
            if other is None:
                grid_opts = [self.check_cell(x, y) for x, y in grid_coords]
            else:
                grid_opts = []
                for n in range(len(grid_coords)):
                    x, y = grid_coords[n]
                    if (x, y) in overlapping_coords.keys():
                        grid_opts.append(self.check_cell(x, y, other=other, other_xy=overlapping_coords[(x, y)]))
                    else:
                        grid_opts.append(self.check_cell(x, y))
            for number in range(1, 10):
                if unlist(grid_opts).count(number) == 1:
                    cell = [n for n, l in enumerate(grid_opts) if number in l][0]
                    x, y = grid_coords[cell]
                    self.sudoku[x, y] = number
                    if other is not None:
                        if (x, y) in overlapping_coords.keys():
                            other.sudoku[overlapping_coords[(x, y)][0], overlapping_coords[(x, y)][1]] = number
                    if plot_progress:
                        if other is not None:
                            self.multiplot(other=other, curr_x=x, curr_y=y, overlapping_coords=overlapping_coords, outfile='prog_plots/{}.png'.format(len(os.listdir('prog_plots'))))
                        else:
                            self.plot(curr_x=x, curr_y=y)

    def check_row(self, x, other=None, overlapping_coords=None, plot_progress=False):
        if np.any(self.sudoku[x, :] == 0):
            if other is None:
                row_opts = [self.check_cell(x, y) for y in range(9)]
            else:
                row_opts = []
                for y in range(9):
                    if (x, y) in overlapping_coords.keys():
                        row_opts.append(self.check_cell(x, y, other=other, other_xy=overlapping_coords[(x, y)]))
                    else:
                        row_opts.append(self.check_cell(x, y))
            for number in range(1, 10):
                if unlist(row_opts).count(number) == 1:
                    y = [n for n, l in enumerate(row_opts) if number in l][0]
                    self.sudoku[x, y] = number
                    if other is not None:
                        if (x, y) in overlapping_coords.keys():
                            other.sudoku[overlapping_coords[(x, y)][0], overlapping_coords[(x, y)][1]] = number
                    if plot_progress:
                        if other is not None:
                            self.multiplot(other=other, curr_x=x, curr_y=y, overlapping_coords=overlapping_coords, outfile='prog_plots/{}.png'.format(len(os.listdir('prog_plots'))))
                        else:
                            self.plot(curr_x=x, curr_y=y)

    def check_col(self, y, other=None, overlapping_coords=None, plot_progress=False):
        if np.any(self.sudoku[:, y] == 0):
            if other is None:
                row_opts = [self.check_cell(x, y) for x in range(9)]
            else:
                row_opts = []
                for x in range(9):
                    if (x, y) in overlapping_coords.keys():
                        row_opts.append(self.check_cell(x, y, other=other, other_xy=overlapping_coords[(x, y)]))
                    else:
                        row_opts.append(self.check_cell(x, y))
            for number in range(1, 10):
                if unlist(row_opts).count(number) == 1:
                    x = [n for n, l in enumerate(row_opts) if number in l][0]
                    self.sudoku[x, y] = number
                    if other is not None:
                        if (x, y) in overlapping_coords.keys():
                            other.sudoku[overlapping_coords[(x, y)][0], overlapping_coords[(x, y)][1]] = number
                    if plot_progress:
                        if other is not None:
                            self.multiplot(other=other, curr_x=x, curr_y=y, overlapping_coords=overlapping_coords,
                                           outfile='prog_plots/{}.png'.format(len(os.listdir('prog_plots'))))
                        else:
                            self.plot(curr_x=x, curr_y=y)

    def check_diag(self, other=None, overlapping_coords=None, plot_progress=False):
        diag1 = [[i, abs(8-i)] for i in range(9)]
        diag2 = [[i, i] for i in range(9)]
        for diag in [diag1, diag2]:
            if other is not None:
                diag_opts = [self.check_cell(x, y) for x, y in diag]
            else:
                diag_opts = []
                for n in range(len(diag)):
                    x, y = diag[n]
                    if (x, y) in overlapping_coords.keys():
                        other_xy = overlapping_coords[(x, y)]
                        diag_opts.append(self.check_cell(x, y, other=other, other_xy=other_xy))
                    else:
                        diag_opts.append(self.check_cell(x, y))
            for number in range(1, 10):
                if unlist(diag_opts).count(number) == 1:
                    cell = [n for n, l in enumerate(diag_opts) if number in l][0]
                    x, y = diag[cell]
                    self.sudoku[x, y] = number
                    if other is not None:
                        if (x, y) in overlapping_coords.keys():
                            other.sudoku[overlapping_coords[(x, y)][0], overlapping_coords[(x, y)][1]] = number
                    if plot_progress:
                        if other is not None:
                            self.multiplot(other=other, curr_x=x, curr_y=y, overlapping_coords=overlapping_coords,
                                           outfile='prog_plots/{}.png'.format(len(os.listdir('prog_plots'))))
                        else:
                            self.plot(curr_x=x, curr_y=y)

    def check_indiv(self, other=None, overlapping_coords=None, plot_progress=False):
        for x in range(9):
            for y in range(9):
                if self.sudoku[x, y] == 0:
                    if plot_progress:
                        if other is not None:
                            self.multiplot(other=other, curr_x=x, curr_y=y, overlapping_coords=overlapping_coords, outfile='prog_plots/{}.png'.format(len(os.listdir('prog_plots'))))
                        else:
                            self.plot(curr_x=x, curr_y=y)
                    if other is not None:
                        if (x, y) in overlapping_coords.keys():
                            cell_opts = self.check_cell(x, y, other=other, other_xy=overlapping_coords[(x, y)])
                            if len(cell_opts) == 1:
                                self.sudoku[x, y] = cell_opts[0]
                                other.sudoku[overlapping_coords[(x, y)][0], overlapping_coords[(x, y)][1]] = cell_opts[0]
                        else:
                            cell_opts = self.check_cell(x, y)
                            if len(cell_opts) == 1:
                                self.sudoku[x, y] = cell_opts[0]
                    else:
                        cell_opts = self.check_cell(x, y)
                        if len(cell_opts) == 1:
                            self.sudoku[x, y] = cell_opts[0]

    def solve(self, max_iter=1000, n_deep=0, max_n_deep=1, plot_progress=False, result_plot=None):
        n = 0
        while (np.any(self.sudoku == 0)) and (n < max_iter):
            start = self.sudoku.copy()
            self.check_indiv(plot_progress=plot_progress)
            _ = [self.check_row(x, plot_progress=plot_progress) for x in range(9)]
            _ = [self.check_col(y, plot_progress=plot_progress) for y in range(9)]
            _ = [self.check_grid(n, plot_progress=plot_progress) for n in range(9)]
            if self.diag:
                self.check_diag(plot_progress=plot_progress)
            if np.all(start == self.sudoku):
                solved = False
                if n_deep < max_n_deep:
                    if 0 in self.sudoku:
                        for x, y in np.array(np.where(self.sudoku == 0)).transpose():
                            if not solved:
                                for opt in self.check_cell(x, y):
                                    self_copy = copy.deepcopy(self)
                                    self_copy.sudoku[x, y] = opt
                                    self_copy.guess_coords += [(x, y)]
                                    solved, s = self_copy.solve(max_iter, n_deep+1, max_n_deep, plot_progress, result_plot)
                if solved:
                    self.sudoku = self_copy.sudoku
                    if result_plot is not None:
                        self.plot(outfile=result_plot)
                    return True, self.sudoku
                else:
                    return False, self.sudoku
            n += 1
        if result_plot is not None:
            self.plot(outfile=result_plot)
        return True, self.sudoku

    def multi_solve(self, other, overlapping_grids, max_iter=1000, n_deep=0, max_n_deep=1, plot_progress=False, result_plot=None):
        if plot_progress:
            if not os.path.isdir('prog_plots'):
                os.makedirs('prog_plots')
        overlapping_coords = {}
        for g1, g2 in overlapping_grids:
            g1_coords = [tuple(x) for x in list(np.array(np.where(self.grid == g1)).transpose())]
            g2_coords = [tuple(x) for x in list(np.array(np.where(other.grid == g2)).transpose())]
            for n in range(len(g1_coords)):
                overlapping_coords[g1_coords[n]] = g2_coords[n]
        inv_overlapping_coords = {v: k for k, v in overlapping_coords.items()}
        n = 0
        while ((np.any(self.sudoku == 0)) or (np.any(other.sudoku == 0))) and (n < max_iter):
            start1 = self.sudoku.copy()
            start2 = other.sudoku.copy()
            self.check_indiv(other=other, overlapping_coords=overlapping_coords, plot_progress=plot_progress)
            other.check_indiv(other=self, overlapping_coords=inv_overlapping_coords, plot_progress=plot_progress)
            _ = [self.check_row(x, other=other, overlapping_coords=overlapping_coords, plot_progress=plot_progress) for x in range(9)]
            _ = [other.check_row(x, other=self, overlapping_coords=inv_overlapping_coords, plot_progress=plot_progress) for x in range(9)]
            _ = [self.check_col(y, other=other, overlapping_coords=overlapping_coords, plot_progress=plot_progress) for y in range(9)]
            _ = [other.check_col(y, other=self, overlapping_coords=inv_overlapping_coords, plot_progress=plot_progress) for y in range(9)]
            _ = [self.check_grid(n, other=other, overlapping_coords=overlapping_coords, plot_progress=plot_progress) for n in range(9)]
            _ = [other.check_grid(n, other=self, overlapping_coords=inv_overlapping_coords, plot_progress=plot_progress) for n in range(9)]
            if self.diag:
                self.check_diag(other=other, overlapping_coords=overlapping_coords, plot_progress=plot_progress)
            if other.diag:
                other.check_diag(other=self, overlapping_coords=inv_overlapping_coords, plot_progress=plot_progress)
            n += 1
            if np.all(start1 == self.sudoku) and np.all(start2 == other.sudoku):
                solved = False
                s1 = self.sudoku
                s2 = other.sudoku
                if n_deep < max_n_deep:
                    if 0 in self.sudoku:
                        for x, y in np.array(np.where(self.sudoku == 0)).transpose():
                            if not solved:
                                if (x, y) in overlapping_coords.keys():
                                    for opt in self.check_cell(x, y, other=other, other_xy=overlapping_coords[(x, y)]):
                                        self_copy = copy.deepcopy(self)
                                        other_copy = copy.deepcopy(other)
                                        self_copy.sudoku[x, y] = opt
                                        self_copy.guess_coords += [(x, y)]
                                        other_copy.sudoku[overlapping_coords[(x, y)][0], overlapping_coords[(x, y)][1]] = opt
                                        other_copy.guess_coords += [(overlapping_coords[(x, y)][0], overlapping_coords[(x, y)][1])]
                                        solved, s1, s2 = self_copy.multi_solve(other_copy, overlapping_grids, max_iter, n_deep + 1, max_n_deep, plot_progress)
                                else:
                                    for opt in self.check_cell(x, y):
                                        self_copy = copy.deepcopy(self)
                                        other_copy = copy.deepcopy(other)
                                        self_copy.sudoku[x, y] = opt
                                        self_copy.guess_coords += [(x, y)]
                                        solved, s1, s2 = self_copy.multi_solve(other_copy, overlapping_grids, max_iter, n_deep + 1, max_n_deep, plot_progress)
                    if (0 in other.sudoku) and (not solved):
                        for x, y in np.array(np.where(other.sudoku == 0)).transpose():
                            if not solved:
                                if (x, y) not in inv_overlapping_coords.keys():
                                    for opt in other.check_cell(x, y):
                                        self_copy = copy.deepcopy(self)
                                        other_copy = copy.deepcopy(other)
                                        other_copy.sudoku[x, y] = opt
                                        other_copy.guess_coords += [(x, y)]
                                        solved, s1, s2 = other_copy.multi_solve(self_copy, [[x[1], x[0]] for x in overlapping_grids], max_iter, n_deep + 1, max_n_deep, plot_progress)
                if solved:
                    other.sudoku = other_copy.sudoku
                    self.sudoku = self_copy.sudoku
                    if result_plot is not None:
                        self.multiplot(other=other, overlapping_coords=overlapping_coords, outfile=result_plot)
                    if plot_progress:
                        self.multiplot(other=other, overlapping_coords=overlapping_coords, outfile='prog_plots/{}.png'.format(len(os.listdir('prog_plots'))))
                    return True, self.sudoku, other.sudoku
                return solved, s1, s2
        if result_plot is not None:
            self.multiplot(other=other, overlapping_coords=overlapping_coords, outfile=result_plot)
        if plot_progress:
            self.multiplot(other=other, overlapping_coords=overlapping_coords, outfile='prog_plots/{}.png'.format(len(os.listdir('prog_plots'))))
        return True, self.sudoku, other.sudoku

    def plot_grid(self, ax, grid_layout=None):
        if grid_layout is None:
            grid_layout = self.grid
        for grid in [x for x in np.unique(grid_layout) if x >= 0]:
            grid_coords = [tuple(x) for x in list(np.array(np.where(grid_layout == grid)).transpose())]
            for x, y in grid_coords:
                if (x - 1, y) not in grid_coords:
                    ax.hlines(abs((x + 1) - grid_layout.shape[0]) + .5, y - .5, y + .5, linewidth=2)
                else:
                    ax.hlines(abs((x + 1) - grid_layout.shape[0]) + .5, y - .5, y + .5, linewidth=.5, colors='grey')
                if (x + 1, y) not in grid_coords:
                    ax.hlines(abs((x + 1) - grid_layout.shape[0]) - .5, y - .5, y + .5, linewidth=2)
                else:
                    ax.hlines(abs((x + 1) - grid_layout.shape[0]) - .5, y - .5, y + .5, linewidth=.5, colors='grey')
                if (x, y - 1) not in grid_coords:
                    ax.vlines(y - .5, abs((x + 1) - grid_layout.shape[0]) - .5, abs((x + 1) - grid_layout.shape[0]) + .5, linewidth=2)
                else:
                    ax.vlines(y - .5, abs((x + 1) - grid_layout.shape[0]) - .5, abs((x + 1) - grid_layout.shape[0]) + .5, linewidth=.5, colors='grey')
                if (x, y + 1) not in grid_coords:
                    ax.vlines(y + .5, abs((x + 1) - grid_layout.shape[0]) - .5, abs((x + 1) - grid_layout.shape[0]) + .5, linewidth=2)
                else:
                    ax.vlines(y + .5, abs((x + 1) - grid_layout.shape[0]) - .5,abs((x + 1) - grid_layout.shape[0]) + .5, linewidth=.5, colors='grey')

    def plot(self, outfile, other=None, overlapping_grids=None, curr_x=None, curr_y=None, col='orange'):
        if other is not None:
            overlapping_coords = {}
            for g1, g2 in overlapping_grids:
                g1_coords = [tuple(x) for x in list(np.array(np.where(self.grid == g1)).transpose())]
                g2_coords = [tuple(x) for x in list(np.array(np.where(other.grid == g2)).transpose())]
                for n in range(len(g1_coords)):
                    overlapping_coords[g1_coords[n]] = g2_coords[n]
            self.multiplot(other=other, overlapping_coords=overlapping_coords, outfile=outfile)
        else:
            fig, ax = plt.subplots(1, 1, figsize=(5, 5))
            ax.set_xlim(-1, 9)
            ax.set_ylim(-1, 9)
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            for x in range(9):
                for y in range(9):
                    if self.sudoku[x, y] != 0:
                        if (x, y) in self.guess_coords:
                            ax.text(y, abs((x + 1) - self.sudoku.shape[0]), str(self.sudoku[x, y]), ha='center', va='center', color='lightgreen')
                        else:
                            ax.text(y, abs((x + 1) - self.sudoku.shape[0]), str(self.sudoku[x, y]), ha='center', va='center')
            self.plot_grid(ax)
            if curr_x is not None:
                rect = Rectangle((curr_y-.5, abs((curr_x + 1) - 9)-.5), 1, 1, edgecolor=col, facecolor='none', linewidth=3)
                ax.add_patch(rect)
            plt.axis('off')
            plt.savefig(outfile, bbox_inches='tight', pad_inches=0, dpi=160)
            plt.close()

    def multiplot(self, other, overlapping_coords, outfile, curr_x=None, curr_y=None, col='orange'):
        sudo2_add_x = list(overlapping_coords.keys())[0][0] - list(overlapping_coords.values())[0][0]
        sudo2_add_y = list(overlapping_coords.keys())[0][1] - list(overlapping_coords.values())[0][1]
        full_sudo = np.zeros((9 + abs(sudo2_add_x), 9 + abs(sudo2_add_y)), dtype=np.uint8)
        full_grid = np.zeros((9 + abs(sudo2_add_x), 9 + abs(sudo2_add_y)), dtype=np.int8) - 1
        if (sudo2_add_y >= 0) and (sudo2_add_x >= 0):
            full_sudo[0:9, 0:9] = self.sudoku
            full_sudo[sudo2_add_x:9 + sudo2_add_x, sudo2_add_y:(9 + sudo2_add_y)] = other.sudoku
            full_grid[0:9, 0:9] = self.grid
            full_grid[sudo2_add_x:9 + sudo2_add_x, sudo2_add_y:(9 + sudo2_add_y)] = (other.grid + 9)
        elif (sudo2_add_y < 0) and (sudo2_add_x < 0):
            sudo2_add_y *= -1
            sudo2_add_x *= -1
            curr_x += sudo2_add_x
            curr_y += sudo2_add_y
            full_sudo[0:9, 0:9] = other.sudoku
            full_sudo[sudo2_add_x:9 + sudo2_add_x, sudo2_add_y:(9 + sudo2_add_y)] = self.sudoku
            full_grid[0:9, 0:9] = other.grid
            full_grid[sudo2_add_x:9 + sudo2_add_x, sudo2_add_y:(9 + sudo2_add_y)] = (self.grid + 9)
        elif (sudo2_add_x >= 0) and (sudo2_add_y < 0):
            sudo2_add_y *= -1
            curr_y += sudo2_add_y
            full_sudo[sudo2_add_x:9 + sudo2_add_x, 0:9] = self.sudoku
            full_sudo[0:9, sudo2_add_y:(9 + sudo2_add_y)] = other.sudoku
            full_grid[sudo2_add_x:9 + sudo2_add_x, 0:9] = self.grid
            full_grid[0:9, sudo2_add_y:(9 + sudo2_add_y)] = (other.grid + 9)
        elif (sudo2_add_x < 0) and (sudo2_add_y >= 0):
            sudo2_add_x *= -1
            curr_x += sudo2_add_x
            full_sudo[sudo2_add_x:9 + sudo2_add_x, 0:9] = other.sudoku
            full_sudo[0:9, sudo2_add_y:(9 + sudo2_add_y)] = self.sudoku
            full_grid[sudo2_add_x:9 + sudo2_add_x, 0:9] = other.grid
            full_grid[0:9, sudo2_add_y:(9 + sudo2_add_y)] = (self.grid + 9)
        fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        ax.set_xlim(-1, full_sudo.shape[1])
        ax.set_ylim(-1, full_sudo.shape[0])
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        for x in range(full_sudo.shape[0]):
            for y in range(full_sudo.shape[1]):
                if full_sudo[x, y] != 0:
                    if (x, y) in self.guess_coords:
                        ax.text(y, abs((x + 1) - full_sudo.shape[0]), str(full_sudo[x, y]), ha='center', va='center', color='lightgreen')
                    else:
                        ax.text(y, abs((x + 1) - full_sudo.shape[0]), str(full_sudo[x, y]), ha='center', va='center')
        self.plot_grid(ax, full_grid)
        if curr_x is not None:
            rect = Rectangle((curr_y - .5, abs((curr_x + 1) - full_sudo.shape[0]) - .5), 1, 1, edgecolor=col, facecolor='none', linewidth=3, zorder=1000)
            ax.add_patch(rect)
        plt.axis('off')
        plt.savefig(outfile, bbox_inches='tight', pad_inches=0, dpi=160)
        plt.close()



if __name__ == '__main__':
    '''
    sudoku_string is all numbers in the sudoku counted from top left to bottom right, missing numbers are 0
    grids should be numbered 0-8 and final gridstrings should also be counted from top left to bottom right
    Here's an example of how these strings would be made:
    
    full figure_like sudoku    sudoku_string per row    grid_string per row
    
     1 2 3 | 7 8 9 | 4 5 6     123789456                000111222
     4 5 6 |       | 7 8 9     456000789                000111222
           |       | 1 2 3     000000123                000111222
     ---------------------
     2 3 1 |       |           231000000                333444555  
           | 2 3 1 | 8 9 7     000231897                333444555
           |       | 2 3 1     000000231                333444555
    ----------------------
     3 1 2 | 9 7 8 |           312978000                333444555
     6   5 |   1   | 9  8      605010908                333444555
       7   | 6   5 |  1        070605010                333444555
    
    So the eventual inputs would be:
    sudoku_string = '123789456456000789000000123231000000000231897000000231312978000605010908070605010'
    grid_string = '000111222000111222000111222333444555333444555333444555666777888666777888666777888'
    This grid_string is the default, so for standard sudokus you dont have to type it out every time.
    
    Enabling the diag argument will force the sudoku to find a solution where the diagonals adhere to the same rules as rows and columns.
    
    For solving connected sudokus: overlapping rids should always be counted from top left to bottom right
    Enabling the plog_progress argument in solve will plot the whole sudoku at each step 
    (orange box: cell being considered, green number: guess) to a directory prog_plots. This will make the solver a LOT slower though.
    
    
    Here are some examples:
    '''
    sudo_in11 = '090802010400030009005000700100000004040000050300000001001000000500070000060704000'
    sudo_gr11 = '000012222011111112033414552003444522033444552333345555666666888777766888777776888'
    sudo_in12 = '000306000000000000000807000809040103000705000201030706000409000000000000000503000'
    sudo_gr12 = '000111222000111222000111222333444555333444555333444555666777888666777888666777888'
    sudo11 = Sudoku(sudo_in11, sudo_gr11)
    sudo12 = Sudoku(sudo_in12, sudo_gr12, diag=True)
    sudo11.plot('sudo1_1.png')
    sudo12.plot('sudo1_2.png')
    sudo11.plot('sudo1.png', other=sudo12, overlapping_grids=[[8, 0]])

    solution1 = sudo11.multi_solve(sudo12, [[8, 0]], max_iter=200, max_n_deep=1, plot_progress=False, result_plot='solution1.png')

    sudo_in21 = '007000900060070030100307005805000709000000000604000103500406000050040000008000000'
    sudo_gr21 = '000111222300112222300112244300115544335555544335566444337766888777766888777666888'
    sudo_in22 = '000000300000010000000209001502000706000000000806000904900305002000040000003000500'
    sudo_gr22 = '000111222000111222000111222333444555333444555333444555666777888666777888666777888'
    sudo21 = Sudoku(sudo_in21, sudo_gr21)
    sudo22 = Sudoku(sudo_in22, sudo_gr22, diag=True)
    sudo21.plot('sudo2_1.png')
    sudo22.plot('sudo2_2.png')
    sudo21.plot('sudo2.png', other=sudo22, overlapping_grids=[[8, 0]])
    solution2 = sudo21.multi_solve(sudo22, [[8, 0]], max_iter=200, max_n_deep=1, plot_progress=False, result_plot='solution2.png')

