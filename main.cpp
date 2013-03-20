/*
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
	name: Alexander Lapin
	e-mail: lapinra@gmail.com
*/

#include <iostream>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <cmath>
#include <algorithm>
#include <utility>

using namespace std;
using std::cin;
using std::cout;
using std::endl;
using std::min;
using std::max;
using std::vector;

typedef vector < double > Point;

const double eps = 1e-8;

double dist(const Point &a, const Point &b)
{
  double result = 0.;
  int n = a.size();

  for (int index = 0; index < n; ++index)
    result += (a[index] - b[index]) * (a[index] - b[index]);

  return result;
}

double norm(const vector < Point > &a)
{
  double result = 0.;
  int n, m;
  n = a.size();
  m = a[0].size();

  for (int first_index = 0; first_index < n; ++first_index)
    for (int second_index = 0; second_index < m; ++second_index)
      result += a[first_index][second_index] * a[first_index][second_index];

  return result;
}

double get_U_ij(const Point &a, const Point &b)
{
  double result = pow(dist(a, b), -6.0) - pow(dist(a, b), -3.0);

  return result;
}

double get_U_tot(const vector < Point > &points)
{
  double result = 0.;
  int n = points.size();

  for (int first_index = 0; first_index < n; ++first_index)
    for (int second_index = first_index + 1; second_index < n; ++second_index)
      result += get_U_ij(points[first_index], points[second_index]);

  return result;
}

double get_U_at(const vector < Point > &points)
{
  int n = points.size();
  return get_U_tot(points) / n;
}

void init_basis_vectors(vector < vector < Point > > &s, int n, int m)
{
  s.assign(n * m, vector < Point > (n, Point(m, 0.)));

  for (int first_index = 0; first_index < n; ++first_index)
    for (int second_index = 0; second_index < m; ++second_index)
    {
      s[first_index * m + second_index][first_index][second_index] = 1.0;
    }
}

vector < Point > normed_residual(vector < Point > &a, vector < Point > &b)
{
  int n = a.size();
  int m = a[0].size();
  vector < Point > c(n, Point(m, 0.));
  for (int first_index = 0; first_index < n; ++first_index)
    for (int second_index = 0; second_index < m; ++second_index)
      c[first_index][second_index] = a[first_index][second_index] - b[first_index][second_index];

  double vect_norm = sqrt(norm(c));
  if (vect_norm > eps)
  {
    for (int first_index = 0; first_index < n; ++first_index)
      for (int second_index = 0; second_index < m; ++second_index)
        c[first_index][second_index] /= vect_norm;
  }

  return c;
}

void add_vector_to_vector_with_coeff(const vector < Point > &add, vector < Point> &res, double lambda)
{
  int n = add.size();
  int m = add[0].size();
  for (int first_index = 0; first_index < n; ++first_index)
    for (int second_index = 0; second_index < m; ++second_index)
      res[first_index][second_index] += lambda * add[first_index][second_index];
}

void get_lambda_segment(const vector < Point > &answer, const vector < Point > &s, double d, double &a, double &c)
{
  double left = 0.;
  double middle = 0.;
  double right = 0.;
  vector < Point > x;
  double y_curr;
  double y_next;

  x = answer;
  y_curr = get_U_tot(x);
  y_next = y_curr;

  double r = (3. - sqrt(5.0)) / 2.;
  r = (1. - r) / r;

  do
  {
    y_curr = y_next;

    left = middle;
    middle = right;
    right += d;

    x = answer;
    add_vector_to_vector_with_coeff(s, x, right);
    y_next = get_U_tot(x);
    
    d = d * r;
  }
  while (y_next < y_curr);

  a = left;
  c = right;
}

double gold_method(const vector < Point > &answer, const vector < Point> &s, double a, double c)
{
  double b, d;
  double r = (3.0 - sqrt(5.0)) / 2.0;
  double y_b, y_d;

  vector < Point > x_d, x_b;

  b = a + r * (c - a);
  x_b = answer;
  add_vector_to_vector_with_coeff(s, x_b, b);
  y_b = get_U_tot(x_b);

  while (abs(c - a) > eps)
  {
    d = c + r * (a - c);
    x_d = answer;
    add_vector_to_vector_with_coeff(s, x_d, d);
    y_d = get_U_tot(x_d);

    if (y_d < y_b)
    {
      a = b;
      b = d;
      y_b = y_d;
    }
    else
    {
      c = a;
      a = d;
    }
  }
  
  return c;
}

double get_best_lambda(const vector < Point > &answer, const vector < Point > &s)
{
  double a, c;
  double lambda_1, lambda_2;
  double delta = 0.000001;
  double best_lambda, best_U;
  double U;
  vector < Point > x;

  get_lambda_segment(answer, s, delta, a, c);
//  cout << a << " " << c << endl;
  lambda_1 = gold_method(answer, s, a, c);

  get_lambda_segment(answer, s, -delta, a, c);
//  cout << a << " " << c << endl;
  lambda_2 = gold_method(answer, s, a, c);

  best_lambda = 0;
  best_U = get_U_tot(answer);

  x = answer;
  add_vector_to_vector_with_coeff(s, x, lambda_1);
  U = get_U_tot(x);

  if (U < best_U)
  {
    best_U = U;
    best_lambda = lambda_1;
  }


  x = answer;
  add_vector_to_vector_with_coeff(s, x, lambda_2);
  U = get_U_tot(x);

  if (U < best_U)
  {
    best_U = U;
    best_lambda = lambda_2;
  }

//  cout << best_lambda << endl;
  return best_lambda;
}

vector < Point > apply_Powell_method(vector < Point > answer, int m)
{
  int n = answer.size();

  int count = 0;
  double U_tot, U_tot_prev;
  double U_at, U_at_prev;
  vector < vector < Point > > s;
  vector < Point > prev_answer;
  init_basis_vectors(s, n, m);

  do
  {
//    U_tot_prev = get_U_tot(answer);
    U_at_prev = get_U_at(answer);

    for (int k = 0; k < 3 * m * n; ++k)
    {
      double lambda;
      prev_answer = answer;
      for (int index = 0; index < n * m; ++index)
      {
        lambda = get_best_lambda(answer, s[index]);

        add_vector_to_vector_with_coeff(s[index], answer, lambda);
        //      cout << "U_tot: " << get_U_tot(answer) << endl;
      }

      for (int index = 0; index < n * m - 1; ++index)
      {
        s[index] = s[index + 1];
      }
      s[n * m - 1] = normed_residual(answer, prev_answer);

      lambda = get_best_lambda(answer, s[n * m - 1]);
      add_vector_to_vector_with_coeff(s[n * m - 1], answer, lambda);

      count++;
      if (count > m * n)
      {
        init_basis_vectors(s, n, m);
        count = 0;
      }
    }

//    U_tot = get_U_tot(answer);
    U_at = get_U_at(answer);

    cout << fabs((U_at - U_at_prev) / U_at) << " " << U_at << endl;
  }
  while (fabs((U_at - U_at_prev) / U_at) > eps);

  return answer;
}

vector < Point > get_gradient(vector < Point > points)
{
  int n = points.size();
  int m = points[0].size();
  vector < Point > result(n, Point(m, 0.));
  double d;

  for (int i = 0; i < n; ++i)
    for (int j = 0; j < m; ++j)
    {
      for (int index = 0; index < n; ++index)
      {
        if (index == i)
          continue;

        d = dist(points[i], points[index]);

//        result[i][j] += 2.0 * (points[i][j] - points[index][j]) * (-6.0 * pow(d, -7.0) - 3.0 * pow(d, -4.0)) / n;
        result[i][j] += -2.0 * (points[i][j] - points[index][j]) * (6.0 - 3.0 * pow(d, 3.0)) * pow(d, -7.0) / n;
//        cout << points[i][j] - points[index][j] << endl;
      }
//      cout << result[i][j] << endl;
    }

  return result;
}

vector < Point > apply_dynamic_method(vector < Point > answer)
{
  int n = answer.size();
  int m = answer[0].size();
  
  double U_at, U_at_prev;

  double dt = 0.085;
  double betta = 0.998;

  vector < Point > v(n, Point(m, 0.));
  vector < Point > grad;

  double c, mod;

  do
  {
    U_at_prev = get_U_at(answer);

    grad = get_gradient(answer);

    if (norm(grad) > 4)
    {
      for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
          grad[i][j] /= 100 * sqrt(norm(grad));
    }

    add_vector_to_vector_with_coeff(grad, v, -dt);
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < m; ++j)
        v[i][j] *= betta;

    add_vector_to_vector_with_coeff(v, answer, dt);

    U_at = get_U_at(answer);
//    cout << fabs((U_at - U_at_prev) / U_at) << " " << U_at << endl;
/*
    for (int i = 0; i < n; ++i)
    {
      for (int j = 0; j < m; ++j)
        cout << answer[i][j] << " ";
      cout << endl;
    }
    cout << endl;
*/
  }
  while (fabs((U_at - U_at_prev) / U_at) > eps);

  return answer;
}

int main()
{
  vector < pair < pair < double, double >, pair < int, int > > > res;
  vector < Point > points;
  FILE *out = fopen("output.txt", "w");

  int n = 5;
  int m = 1;
  for (int n = 3; n < 20 + 1; ++n)
    for (int m = 1; m < 6 + 1; ++m)
    {
      points.assign(n, Point(m, 0.));
      for (int first_index = 0; first_index < n; ++first_index)
        for (int second_index = 0; second_index < m; ++second_index)
        {
          points[first_index][second_index] = -3.0 + 6.0 * (1.0 * rand() / RAND_MAX);
//          cout << points[first_index][second_index] << endl;
        }

      cout << "U initial: " << get_U_tot(points) << endl;

//      points = apply_Powell_method(points, m);
      cout << "Starting dynamic:" << endl;
      points = apply_dynamic_method(points);
      cout << "Starting Powell:" << endl;
      points = apply_Powell_method(points, m);
      cout << "Starting dynamic:" << endl;
      points = apply_dynamic_method(points);
      cout << "Starting Powell:" << endl;
      points = apply_Powell_method(points, m);

      for (int i = 0; i < n; ++i)
      {
        for (int j = 0; j < m; ++j)
          cout << points[i][j] << " ";
        cout << endl;
      }
      cout << "n = " << n << " m = " << m << "; U_tot: " << get_U_tot(points) << "; U_at: " << get_U_at(points) << endl;

      res.push_back(make_pair(make_pair(get_U_at(points), get_U_tot(points)), make_pair(n, m)));
    }

  sort(res.begin(), res.end());
  for (int i = 0; i < res.size(); ++i)
    fprintf(out, "%d %d %.4lf %4lf\n", res[i].second.first, res[i].second.second, res[i].first.second, res[i].first.first);

  fclose(out);

  return 0;
}

