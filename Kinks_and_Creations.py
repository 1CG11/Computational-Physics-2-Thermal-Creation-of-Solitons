"""
Defines functions to identify kink anti-kink pairs and creations 
as described in sections 4 and 5:

zero_crossings
zeros_and_wide_gaps
kink_in_block
anti_kink_in_block
pairs
smooth
creations
creation_rates

Defines the variables:
w_kink, h_kink, d_kink, d_kink_frame, buff_frame
"""
from Discretisation import np, N, dx, dt, frame_space
from Initial_Conditions import heat_bath
                         
# kink variables
w_kink = 20                 # minimum width of kink (nodes)
h_kink = 0.5                # minimum height of a kink 
d_kink = 50                 # minimum duration of kink (time-steps)

# can also be expressed in frames
d_kink_frame = d_kink // frame_space

# a buffer period required to use smooth properly
buff_frame = d_kink_frame * (d_kink_frame - 1)//2


def zero_crossings( f ):
    """
    Identifies the zero-crossings of 'f'
    Returns the indices of the array 'f' where the array value 
    differs in sign to its predecessor
    
    Assumes array is non-zero
    """
    return np.where( f * np.roll(f,1) < 0)[0]


def zeros_and_wide_gaps( f ):
    """
    Returns number of zero-crossings and wide gaps
    as described in section 4.2
    using 'zero_crossings''
    """
    zeros = zero_crossings( f )
    
    #   if no zero-crossings, return no wide gaps
    z = len(zeros)
    if z == 0:
        return 0, 0
    
    #   else split into blocks
    #   search for a block that is a wide gap
    gap_count = 0
    
  
    #   define boundary section
    block = np.append( np.arange(zeros[-1], N), np.arange(0, zeros[0]))
    
    # if passes width requirement, increase count
    if len(block) >= w_kink:
        gap_count += 1
            
    #   every other section
    for i in range( 1, z ):
        
        #   define section
        block = np.arange(zeros[i-1], zeros[i])
        
        # if passes width requirement, increase count
        if len(block) >= w_kink:
            gap_count += 1

    return z, gap_count
    
    
def kink_in_block( block, f ):
    """
    Determines if there is a kink present in the block of 'f';
    if it passes the width and height conditions.
    
    Block: indices of 'f' where the sign of 'f' is constant
    f : the current field configuration
    """
    high_count = 0
    
    # for each node
    for i in block :
        
        # if passes height requirement, then increase counter
        if f[i] > h_kink:
            high_count += 1
            
            # if passes width requirement, then there is a kink
            if high_count >= w_kink:
                return True
            
        # if fails height requirement, then reset counter
        else:
            
            # reset counter
            high_count = 0
            
    # if at end no kink was detected, then there is no kink   
    return False
        

def anti_kink_in_block( block, f ):
    """
    Determines if there is an anti-kink present in the block of 'f';
    if it passes the width and height conditions.
    
    Block: indices of 'f' where the sign of 'f' is constant
    f : the current field configuration
    """
    deep_count = 0
    
    # for each node
    for i in block :
        
        # if passes height requirement, incerase counter
        if f[i] < - h_kink:
            deep_count += 1
            
            # if passes width requirement, then there is an anti - kink 
            if deep_count >= w_kink:
                return True
            
        # if fails height requirement, reset counter
        else:
            deep_count = 0
            
    # if at end no anti - kink was detected, then there is no anti - kink
    return False


def pairs( f ):
    """
    Calculates the pair number 'n'
    Using the procedure described in 4.3 
    using 'zero_crossings', 'kink_in_block' and 'anti_kink_in_block'

    """
    # identify all zero-crossings
    zeros = zero_crossings( f )
    
    #   if no zero-crossings, return no pairs
    z = len(zeros)
    if z == 0:
        return 0
    
    #   else split into blocks
    #   search for a block containing a kink
    
    # set counters to zero
    kink_count = 0
    anti_kink_count = 0
    
    # we are not yet looking for anything in particular
    kink_search = False
    anti_kink_search = False

    for i in range( z ):
        
        # define block
        
        if i == 0:
            # due to periodic boundary conitions,
            # boundary block is defined differently 
            block = np.append( np.arange(zeros[-1], N), \
                              np.arange(0, zeros[0]))
        else:
            block = np.arange(zeros[i-1], zeros[i])
            
            
        if kink_search == True:
            #   search for kink in block
            
            if kink_in_block( block, f ) == True:
                #   if there is, count it, and now search for anti-kink
                kink_count += 1
                kink_search = False
                anti_kink_search = True
                
        if anti_kink_search == True:
            #   search for anti-kink in block
            
            if anti_kink_in_block( block, f ) == True:
                #   if there is, count it, and now search for kink
                anti_kink_count += 1
                kink_search = True
                anti_kink_search = False
                
        else:
            # search for both
            # only happens until first detection
            
            if kink_in_block( block, f ) == True:
                #   if there is, count it, and now search for anti-kink
                kink_count += 1
                anti_kink_search = True
                
            if anti_kink_in_block( block, f ) == True:
                #   if there is, count it, and now search for kink
                anti_kink_count += 1
                kink_search = True
 
    #   return the minimum; the number of pairs
    return min( kink_count, anti_kink_count )

#   Kink-count Smoothing
def smooth(array, tmax_frame):
    """
    Removes fluctuations of duration less than 'd_min' from an array
    according to the procedure described in section 5.1
    """
    #   for each duration 'd' up to the minimum acceptable
    for d in range(1, d_kink_frame ):
        #   remove, left to right, all fluctuations of duration 'n'
        
        #   for each node before the respective buffer
        for k in range(1, len(array) - d * (d + 1)//2 ):
            #   if a value doesn't match its predecessor
            if array[k] != array[k - 1]:
                #   set it to the value 'n' nodes in the future
                array[k] = array[k + d]
                
    return array[: tmax_frame ]
    
    
def creations(array, tmax_frame):
    """
    Calculates the number of creations that occurred in this time-frame
    given the array of pair numbers over (tmax_frame + buff_frame) frames
    using the procedure described in section 5.1
    using 'smooth'
    """
    s_array = smooth(array, tmax_frame)
    #   count the upward steps
    diff = [max(y - x, 0) for x, y in zip(s_array[:-1], s_array[1:])]
    amount = sum( diff )
    return amount


def creation_rates( pair_array, tmax_frame ):
    """
    Calculates the creation time 'tau' and creation rate 'Gamma'
    given the array of pair numbers over (tmax_frame + buff_frame) frames
    using the procedure described in section 5.1
    using 'creations'
    """
    # we measure these quantities in proper time
    # 1 frame is 'frame_space' timesteps
    # 1 timestep is 'dt' time units
    # therefore one frame is 'frame_space' * dt time units
    tmax_units = frame_space * dt * tmax_frame
    creation_amount = creations(pair_array, tmax_frame)
    
    tau =  tmax_units / creation_amount
    Gamma = creation_amount / tmax_units 
    return Gamma, tau