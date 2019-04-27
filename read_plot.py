import math

def read_plot(filename) :
    try:
        file = open( filename, 'r' )
        lines = file.readlines()
        x = []
        y = []
        z = []
        v = []
        for line in lines :
            if line != '\n' :
                words = line.split()
                ix, iy, iz, iv = [float(s) for s in words]
                x.append( ix )
                y.append( iy )
                z.append( iz )
                v.append( iv )
        return x,y,z,v
    except:
        print 'Error in read_plot. Returning None'
        return None
