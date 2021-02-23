from pylab import *

def angle_between(v1,v2):
    angle = np.arctan2(v2[1],v2[0]) - np.arctan2(v1[1],v1[0])
    return(angle)

def accel_for_force(com,f,r,mass,moi):    
    """Returns the linear and angular accelerations that result from
    the application of a force f, at point r, to an object whose
    centre of mass is com. 

    com  -- centre of mass of accelerated object
    f    -- force vector
    r    -- point at which force is applied
    mass -- mass of accelerated object
    moi  -- moment of inertia of accelerated object

    """
    
    #### linear acceleration
    # f = ma; a = f/m
    ax,ay = f / mass

    #### torque
    # c: vector from centre to application point
    c = r-com[:2]

    if np.linalg.norm(c) == 0.0 :
        # force is applied at the centre of mass. There is no torque.
        at = 0.0
    else:
        ## angle between c and f
        theta = angle_between(c,f) #np.arccos( cos_theta )

        ## torque = r |f| sin theta, but r here is the magnitude of r
        t = np.linalg.norm(r) * np.linalg.norm(f) * np.sin(theta)
        at = t / moi
    
    return ax,ay,at


if __name__ == '__main__' :
    ## centre of mass
    mass = 0.05 # units
    ## moment of inertia
    moi = 0.05 # units 
    p = x,y,a = np.array([0.0,0.0,np.pi/2])
    v = np.array([0.0,0.0])
    r = 1.0 # half-length of rod

    ## force vector
    f = np.array([1.0,0.0])

    ## point where force is applied


    r = np.array([r*np.cos(a),r*np.sin(a)]) # tip of rod
    # r = np.array([0.0*r*np.cos(a),0.0*r*np.sin(a)]) # centre of rod

    #r = np.array([-0.1*r*np.cos(a),-0.1*r*np.sin(a)]) # to right of centre
    #r = np.array([+0.1*r*np.cos(a),+0.1*r*np.sin(a)]) # to left of centre of rod

    dvx,dvy,dva = accel_for_force(p,f,r,mass,moi)

    ## instantaneous accel, changes velocity once/isntantly
    vx = dvx
    vy = dvy
    va = dva

    print(vx,vy,va)

    DT = 0.01

    for it in range(25):
        x,y,a = p
        plot(x,y,'.',alpha=0.5,color='0.0')

        r = 1.0
        line_x = [x+cos(a)*-r,x+cos(a)*r]
        line_y = [y+sin(a)*-r,y+sin(a)*r]
        plot(line_x,line_y,color='0.0',alpha=0.2)
        p += np.array([vx,vy,va]) * DT

    show()
