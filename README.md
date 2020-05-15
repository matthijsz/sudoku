# Sudokusolver

This is my Python sudoku solver, it can do standard Sudokus, Sudokus with a different grid layout, and 2 connected sudokus.

# Usage

### Formatting the input

sudoku_string is all numbers in the sudoku counted from top left to bottom right, missing numbers are 0
grids should be numbered 0-8 and final gridstrings should also be counted from top left to bottom right
Here's an example of how these strings would be made:
```
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
grid_string = '000111222300112222300112244300115544335555544335566444337766888777766888777666888'
```
This standard grid is the default for grid_string, so if you are only solving a standard sudoku you don't have to specify this.

### Diagonal

Enabling the `diag` argument will force the sudoku solver to find a solution where the diagonals adhere to the same rules as rows and columns (contain all numbers 1-9).

### Single sudokus

For solving single sudokus first create the Sudoku object with your input:
```python
from sudokusolver import Sudoku
sudo_in = '000306000000000000000807000809040103000705000201030706000409000000000000000503000'
sudo_gr = '000111222000111222000111222333444555333444555333444555666777888666777888666777888'
my_sudo = Sudoku(sudo_in, sudo_gr)
```
Optionally plot the sudoku to verify your input
```python
my_sudo.plot('sudo1_1.png')
```
Then solve the sudoku
```python
my_solution = my_sudo.solve(max_iter=200, max_n_deep=1, plot_progress=False, result_plot='my_sudoku_solution')
```
The solve function will return a numpy array containing the solution.

### Multiple sudokus

For solving 2 connected sudokus first create a Sudoku object for both with input as described above.

Then run `multisolve` for one of them, and specify `overlapping_grids`. This should be a nested list with length equal to the number of overlapping grids, each element is a list with `[grid_sudoku1, grid_sudoku2]`. 
If there are multiple overlapping grids add them to the list from top-left to bottom-right
For example: [[8, 0]] if only the bottom-right grid in the first sudoku and top-left grid in the second sudoku overlap.

To solve multiple sudokus first create separate Sudoku objects for both:
```python
from sudokusolver import Sudoku
sudo_in11 = '090802010400030009005000700100000004040000050300000001001000000500070000060704000'
sudo_gr11 = '000012222011111112033414552003444522033444552333345555666666888777766888777776888'
sudo_in12 = '000306000000000000000807000809040103000705000201030706000409000000000000000503000'
sudo_gr12 = '000111222000111222000111222333444555333444555333444555666777888666777888666777888'
sudo11 = Sudoku(sudo_in11, sudo_gr11)
sudo12 = Sudoku(sudo_in12, sudo_gr12, diag=True)
```
Optionally plot the sudokus to verify your input
```python
sudo11.plot('sudo1_1.png')
sudo12.plot('sudo1_2.png')
sudo11.plot('sudo1.png', other=sudo12, overlapping_grids=[[8, 0]])
```
Then solve the multi-sudoku
```python
my_solutions = sudo11.multi_solve(sudo12, [[8, 0]], max_iter=200, max_n_deep=1, plot_progress=False, result_plot='solution1.png')
```
multi_solve will return a list of two numpy arrays containing the solutions for the left and right sudoku

### Other arguments
`max_iter` _(int)_ defines the maximum number of iterations the solver will run. This is generally not necessary as it will stop automatically if it can't find a solution.

`max_n_deep` _(int)_ defines how many numbers the solver is allowed to guess. If this is set to 1, and the solver cannot find a solution with the standard strategy it wil iterate over all possible numbers, guessing them, until it finds a solution.
If this is set to two the solver will guess the a number, try to solve, then guess a second number, if it can't solve it, the second guess is changed, untill there are no second guesses left, etc.

`plot_progress` _(bool)_  will create a directory called `prog_plots` and save a full plot at each step (named 0.png, 1.png... etc) in that directory. The orange box is the cell under consideration, and green numbers are guessed numbers.

`result_plot` _(str)_ will plot the solution to this filename.

`n_deep` _(int)_ is only used internally to keep track of how many guesses have been made as it uses the solve or multisolve function recursively.

### Example output
Here's an example gif of a solution. <br/> 

![](out.gif)