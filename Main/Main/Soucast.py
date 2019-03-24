import numpy
class Soucast(object):
    """description of class"""
    T0 = 20

    T = numpy.asarray([T0,20])
    E = numpy.asarray([0,0])
    e = 0
    alfa = numpy.asarray([0,0])

    def lin_interpolace(X,v_X,v_Y):
        #serazeni v_X a v_Y podle v_X
        arr1inds = v_X.argsort()
        sorted_arr1 = v_X[arr1inds[::1]] # -1 by to seradila sestupne, muze se nekdy hodit.
        sorted_arr2 = v_Y[arr1inds]
        return numpy.interp(X,sorted_arr1,sorted_arr2)