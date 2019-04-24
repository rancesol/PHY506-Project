import math

def read_plot(filename) :
    try:
        file = open( filename, 'r' )
        lines = file.readlines()
        x = []
        y = []
        z = []
        for line in lines :
            if line != '\n' :
                words = line.split()
                ix, iy, iz = [float(s) for s in words]
                x.append( ix )
                y.append( iy )
                z.append( iz )
        return x,y,z
    except:
        print 'Error in read_plot. Returning None'
        return None
