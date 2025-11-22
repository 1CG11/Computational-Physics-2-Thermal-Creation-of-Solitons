"""
Produces Data for Table in section 2.1

Comparison of Speeds of 3 implementations of the finite difference method.
"""
from Discretisation import np, C_1,C_2, C_3, N, next_timestep
from Initial_Conditions import heat_bath
import time


T = 1.0
tmax = 10**5


def Explicit_Loop( f_old , f ):
    """
    Given the previous and current field configurations, 'f_old' and 'f',
    updates the field configurations to the next timestep
    according to equation (13)
    
    In particular uses the Explicit Loop
    
    Returns
    -------
    f :     current field configuration
    f_new : next field configuration 
    """
    f_new = - f_old + C_3 * f ** 3
    for i in range ( N ) :
        f_new [ i - 1 ] += C_1 * f [ i - 1 ] + C_2 * ( f [ i - 2 ] + f [ i ] )
    return f , f_new


M = np . zeros ( N * N ) . reshape (( N , N ) )
for i in range( N ):
    M [i , i - 1 ] = C_2
    M [i , i ] = C_1
    M [ i -1 , i ] = C_2
    
def Matrix( f_old , f ):
    """
    Given the previous and current field configurations, 'f_old' and 'f',
    updates the field configurations to the next timestep
    according to equation (13)
    
    In particular uses the Matrix Method equation (14)
    
    Returns
    -------
    f :     current field configuration
    f_new : next field configuration 
    """
    return f , - f_old + np . matmul (M , f ) + C_3 * f * f * f


# test of explicit loop
f_old, f = heat_bath(T)
t0 = time.time()
for _ in range(tmax):
    f_old, f = Explicit_Loop( f_old , f )
t1 = time.time()
t = t1-t0

print('Explicit Loop ' + 
      'Time Taken: '+str(t) + 
      ' Updates per second: '+ str(tmax/t) )



# test of Matrix Method
f_old, f = heat_bath(T)
t0 = time.time()
for _ in range(tmax):
    f_old, f = Matrix( f_old , f )  
t1 = time.time()
t = t1-t0
print('Matrix ' + 
      'Time Taken: '+str(t) + 
      ' Updates per second: '+ str(tmax/t) )


# test of Rolling Array

f_old, f = heat_bath(T)
t0 = time.time()
for _ in range(tmax):
    f_old, f = next_timestep( f_old , f )  
t1 = time.time()
t = t1-t0
print('Rolling Array ' + 
      'Time Taken: '+str(t) + 
      ' Updates per second: '+ str(tmax/t) )