import math

def read_plot(filename) :
    try:
        file = open( filename, 'r' )
        lines = file.readlines()
        x = []
        y = []
        for line in lines :
            if line != '\n' :
                words = line.split()
                ix, iy = [float(s) for s in words]
                x.append( ix )
                y.append( iy )
        return x,y
    except:
        print 'Error in read_plot. Returning None'
        return None
