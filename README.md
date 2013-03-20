Overview
--------

  This program provides Dynamic (gradient) and Pawlle's (not gradient) methods 
  for multi-variable function optimization

  Consider next problem:
  There are `n` points in `m`-dimension space. 
  The potential energy of pair of poits `(p_i, p_j)` is:
    U_ij = (1/r_ij)^12 - (1/r_ij)^6
  Total energy:
    U_tot = Sum U_ij, i < j
  Average energy for point:
    U_at = 1/n U_tot

  We want to minimaze `U_at`.
  Considering `n m` vector as an argument for function, we are looking for best coordinates of points.

  Actually, it finds only local minimums. For better results both methods are used.

  `output.txt` is a results of applying both methods
  `output1.txt` is a results of Pawlle's method
  `output2.txt` is a results of Dynamic method

Some information about methods
------
  Pawlle's method:
  s = {e1, e2, ..., e_q} -- set of basis vectors (for example, q = n m in this task)
  answer -- current answer
  prev_answer -- answer on the previous step

  do {
  prev_answer = answer
  for i = 1 to q
    lambda = argmin_lambda(answer + lambda * s_i)
    answer = answer + lambda * s_i

  for i = 1 to q-1
    s_i = s_{i+1}

  s_q = (answer - prev_answer)
  s_q /= norm(s_q)
  lambda = argmin_lambda(answer + lambda * s_i)
  answer = answer + lambda * s_i

  step++
  if (step > q)
    s = {e1, e2, ..., e_q}
  } while (abs(answer - prev_answer) / answer > eps)

  Moreover it uses some tricks with loops for better converge.


Usage:
------

1. `make`
2. output.txt - for results

License
-------
    Copyright (C) 2013  Alexander Lapin

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
	
Contacts:
---------
Alexander Lapin, lapinra@gmail.com
