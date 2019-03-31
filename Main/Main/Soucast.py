import numpy
class Soucast(object):
    """description of class"""
    T0 = 20

    T = numpy.asarray([T0,0])
    T_Ezk = numpy.asarray([0,0])
    E_zk = numpy.asarray([0.01,0.01])
    e = 0
    T_azk = numpy.asarray([0,0])
    alfa_zk = numpy.asarray([0,0])

    def lin_interpolace(X,v_X,v_Y):
        #serazeni v_X a v_Y podle v_X
        arr1inds = v_X.argsort()
        sorted_arr1 = v_X[arr1inds[::1]] # -1 by to seradila sestupne, muze se nekdy hodit.
        sorted_arr2 = v_Y[arr1inds]
        return numpy.interp(X,sorted_arr1,sorted_arr2)

    def double_interp(T,Q,x,q,t):
        X = numpy.zeros(len(t))
        for temp in range(len(t)):
            #if (Q < min(q[temp])):
            #    lin_k = (numpy.polyfit(x[temp],q[temp],1))
            #    X[temp] = (Q - lin_k[1]) / lin_k[0]                
            #else:
            X[temp] = Soucast.lin_interpolace(Q,q[temp],x[temp])
        X_fin = Soucast.lin_interpolace(T,t,X)
        return X_fin

    def calc_mat_parameter(T,T_zk,Y_zk):
        Y = numpy.zeros([len(T)])
        for t_i,t in enumerate(T):
            Y[t_i] = Soucast.lin_interpolace(t,T_zk,Y_zk)
        return Y