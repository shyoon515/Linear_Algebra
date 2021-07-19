# Lagrange Interpolation Formula

import matplotlib.pyplot as plt
from matplotlib.image import imread
import numpy as np

# class Polynomial: when provided with list of solutions of f(x)=0 and biggest coefficient, makes polynomial f(x). When solution list has no element, the second parameter refers to constant term.
class Polynomial:
    def __init__(self, list_of_x_points=[], biggest_coefficient=1):
        self.coefficients = Polynomial.get_coefficients(list_of_x_points, biggest_coefficient)
        self.dim = len(list_of_x_points)
        self.list_of_x_points = list_of_x_points
        self.biggest_coefficient = biggest_coefficient

    def get_coefficients(list_of_x_points, biggest_coefficient):
        dim = len(list_of_x_points)
        coefficients = [0]*(dim+1)
        for cur_dim in range(dim+1):
            coefficients[dim-cur_dim] = Polynomial.get_coefficient(list_of_x_points, dim, cur_dim)
        for i in range(len(coefficients)):
            coefficients[i] *= biggest_coefficient
        return coefficients

    def get_coefficient(list_of_x_points, dim, cur_dim):
        if cur_dim == dim:
            return 1
        else:
            indices_list = combination(dim, dim-cur_dim)
            summation = 0
            for indices in indices_list:
                multiple = 1
                for index in indices:
                    multiple *= list_of_x_points[index]
                summation += multiple
            
            if dim%2 == 0:    # process of determining signs for each coefficients.
                if cur_dim%2 == 1:
                    summation *= -1
            else:
                if cur_dim%2 == 0:
                    summation *= -1
            return summation

    def __str__(self):
        string_list = ["f(x) =",]
        dim = len(self.coefficients)-1
        for i, cof in enumerate(self.coefficients):
            if cof == 0:
                continue
            elif cof > 0:
                if i == dim:
                    string = " + %g" % cof
                elif cof == 1:
                    string = " + (x^%d)" % (dim-i)
                else:
                    string = " + %g(x^%d)" % (cof, dim-i)
            elif cof < 0:
                if i == dim:
                    string = " - %g" % abs(cof)
                elif cof == -1:
                    string = " - (x^%d)" % (dim-i)
                else:
                    string = " - %g(x^%d)" % (abs(cof), dim-i) 
            string_list.append(string)
        if (self.coefficients[0] > 0) and (self.coefficients[0] != 1) and (dim != 0):
            string_list[1] = " %g(x^%d)" % (self.coefficients[0], dim)
        elif (self.coefficients[0] > 0) and (self.coefficients[0] == 1) and (dim != 0):
            string_list[1] = " (x^%d)" % (dim)
        return "".join(string_list) if len(string_list) != 1 else "f(x) = 0"
    
    def __add__(self, other):
        new_dim = self.dim if self.dim >= other.dim else other.dim
        smaller_dim = self.dim if self.dim <= other.dim else other.dim
        list_of_x_points = [0]*new_dim
        temp_pol = Polynomial(list_of_x_points, 1)
        temp_pol.dim = new_dim
        
        for i in range(smaller_dim+1):
            temp_pol.coefficients[new_dim-i] = self.coefficients[self.dim-i] + other.coefficients[other.dim-i]
        if new_dim != smaller_dim:
            for i in range(new_dim-smaller_dim):
                temp_pol.coefficients[i] = self.coefficients[i] if self.dim > other.dim else other.coefficients[i]
        return temp_pol

    def __sub__(self, other):
        new_dim = self.dim if self.dim >= other.dim else other.dim
        smaller_dim = self.dim if self.dim <= other.dim else other.dim
        list_of_x_points = [0]*new_dim
        temp_pol = Polynomial(list_of_x_points, 1)
        temp_pol.dim = new_dim
        
        for i in range(smaller_dim+1):
            temp_pol.coefficients[new_dim-i] = self.coefficients[self.dim-i] - other.coefficients[other.dim-i]
        if new_dim != smaller_dim:
            for i in range(new_dim-smaller_dim):
                temp_pol.coefficients[i] = self.coefficients[i] if self.dim > other.dim else other.coefficients[i]*(-1)
        return temp_pol


# length of list and the number of picks are the input. returns list of lists, which refers to the indices of combination.
def combination(list_len, pick, recur=0, temp_list=[], index_log=[]):
    if pick == (recur):
        return temp_list
    else:
        re_temp_list = temp_list.copy()
        re_index_log = index_log.copy()
        temp_index_list = []
        for i in range(list_len-recur-sum(index_log)):
            re_index_log.append(i)
            re_temp_list.append(i+recur+sum(index_log))
            occasion = combination(list_len, pick, recur=recur+1, temp_list = re_temp_list, index_log=re_index_log)
            re_temp_list = temp_list.copy()
            re_index_log = index_log.copy()
            if pick == (recur+1):
                temp_index_list.append(occasion)
            else:
                for occ in occasion:
                    temp_index_list.append(occ)
        return temp_index_list
    return index_list


def lagrange_polynomial(list_of_x_points=[]):
    pol_list = []
    for x in list_of_x_points:
        new_x_list = list_of_x_points.copy()
        new_x_list.remove(x)
        big_cof = 1
        for x_2 in new_x_list:
            big_cof *= x-x_2
        pol_list.append(Polynomial(new_x_list, (1/big_cof)))
    return pol_list

def lagrange_interpolation_formula(list_of_points=[]):
    list_of_x_points = []
    list_of_y_points = []
    for x, y in list_of_points:
        list_of_x_points.append(x)
        list_of_y_points.append(y)
    pol_list = lagrange_polynomial(list_of_x_points)
    cp_pol_list = pol_list.copy()

    for i in range(len(cp_pol_list)):
        cp_pol_list[i].biggest_coefficient *= list_of_y_points[i]
        for j in range(len(cp_pol_list[i].coefficients)):
            cp_pol_list[i].coefficients[j] *= list_of_y_points[i]
    lag_pol = cp_pol_list.pop(0)
    for pol in cp_pol_list:
        lag_pol += pol
    return lag_pol


def f(x, polynomial):
    result = 0
    for i in range(polynomial.dim+1):
        result += polynomial.coefficients[i]*x**(polynomial.dim-i)
    return result


def plot_lif(list_of_points):
    list_of_x_points = []
    list_of_y_points = []
    for x, y in list_of_points:
        list_of_x_points.append(x)
        list_of_y_points.append(y)
    
    minimum = int(min(list_of_x_points))
    maximum = int(max(list_of_x_points))
    gap = maximum - minimum
    
    
    if gap >= 6:
        start = minimum-(int(gap/6))
        end = maximum+(int(gap/6))
    else:
        start = minimum-1
        end = maximum+1
    
    lag = lagrange_interpolation_formula(list_of_points)
    lag_show = lag.__str__()
    print(lag)
    x_val = np.linspace(start, end, gap*50)
    y_val = list(map(lambda x : f(x, lag), x_val[:]))
    
    plt.title("Lagrange Interpolation Formula")
    plt.plot(x_val, y_val, label='Interpolation Formula')    # use this when you don't want to see the polynomial.
    #plt.plot(x_val, y_val, label='Interpolation Formula: \n    %s' % lag_show)    #use this when you want to see the polynomial.
    plt.scatter(list_of_x_points, list_of_y_points, c='r', label='Original Input Point')
    plt.grid()
    plt.legend()
    plt.savefig('/workspace/Linear_Algebra/test.png')


plot_lif([(-2,3),(-1,-6),(1,0)])
