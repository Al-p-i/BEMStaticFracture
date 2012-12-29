import math
class Data:
    """data"""
    #to be input from input.csv
    G = 0.0
    v = 0.0
    N = 0
    pressure = []
    fromPoint = 0.0
    toPoint = 0.0

    #to be counted
    length = 0.0
    elemLength = 0.0
    multipleCoef = 0.0
    elementX = []
    arrayD = []
    
class CorelationTableIntensityCoefficient:
    relativeOffset = []
    relativeIntensityCoefficient = []


class TableSigma:
    def __init__(self):
        self.arraySigma = []
        self.X = []
        self.Y = []
        self.charInfo = "NULL"



def readInput(fileName):
    """read input csv"""
    import csv
    data = Data()
    with open(fileName, 'rb') as csvfile:
        inputTable = [row for row in csv.reader(csvfile, delimiter=';')]
        data.G = float(inputTable[0][1])
        data.v = float(inputTable[1][1])
        data.N = int(inputTable[2][1])
        data.pressure = inputTable[3]
        data.pressure = data.pressure[1:]
        for i in range(0, len(data.pressure)):
            data.pressure[i] = float(data.pressure[i])
        
        data.fromPoint = float(inputTable[4][1])
        data.toPoint = float(inputTable[5][1])
        return data


def countInfluenceCoefficient(multipleCoef, elemLength, elementX1, elementX2):
    """count Aij"""
    return multipleCoef * (elemLength / 2.0) / ((elementX1 - elementX2)**2 - (elemLength / 2.0)**2)


def countSpecianInfluenceCoefficient(multipleCoef, elemLength, elementX1, elementX2):
    """count special Aij"""
    s = elemLength / 2.0 + abs(elementX2 - elementX1)
    return  multipleCoef/2.0 * (math.sqrt(2.0) / (s - elemLength) + 1.0 / (2.0 * math.sqrt(elemLength * s / 2.0)) * math.log(abs(math.sqrt(s) - math.sqrt(elemLength)) / (math.sqrt(s) + math.sqrt(elemLength))))


def fillSpecialInfluenceCoefficientMatrix(multipleCoef, elemLength, N, elementX):
    specialCoefficientMatrix = []
    for i in range(0, N):
        specialCoefficientMatrix.append([])
        for j in range(0, N):
            if(j == 0 or j == N - 1):
                specialCoefficientMatrix[i].insert(j, countSpecianInfluenceCoefficient(multipleCoef, elemLength, elementX[i], elementX[j]))
            else:
                specialCoefficientMatrix[i].insert(j, countInfluenceCoefficient(multipleCoef, elemLength, elementX[i], elementX[j]))
    return specialCoefficientMatrix


def fillElementX(N, elemLength, fromPoint):
    elementX = []
    for i in range(0, N):
        elementX.insert(i, fromPoint + i * elemLength + elemLength / 2.0)
    return elementX


def writeDataForSlau(fileName, N, specialCoefficientMatrix, pressure):
    """prepare system.csv"""
    import csv
    outfile = open(fileName, 'wb')
    writer = csv.writer(outfile)

    # write data.N
    nArray = []
    nArray.insert(0, N)

    writer.writerow(nArray)

    # write specialCoefficientMatrix
    for row in specialCoefficientMatrix:
        writer.writerow(row)

    # write negPressure
    negPressure = []
    for i in range(0, len(pressure)):
        negPressure.insert(i, -pressure[i])
        
    writer.writerow(negPressure)

    #initialApproximation = 0.0
    initialApproximation = []
    for i in range(0, len(pressure)):
        initialApproximation.insert(i, 0.0)

    writer.writerow(initialApproximation)


def execSlau(fileName):
    ##TODO wait for slau.exe finish writing out.csv
    cmd = 'slau.exe'
    import subprocess
    #PIPE = subprocess.PIPE
    p = subprocess.Popen(cmd, shell = True)
    p.wait()
    print '\nslau exec finished'


def readOut(fileName):
    """read input csv"""
    import csv
    with open(fileName, 'rb') as csvfile:
        #arrayD = [row for row in csv.reader(csvfile, delimiter=';')]
        #for i in range(0, len(arrayD)):
        #    arrayD.append()
        arrayD = []
        for row in csv.reader(csvfile):
            for value in row:
                arrayD.append(value)

        for i in range(0, len(arrayD)):
            arrayD[i] = float(arrayD[i])

        return arrayD


def countDeformationEnergy(multipleCoef, elemLength, N, elementX, pressure):
    #1 count special coefficient matrix
    specialCoefficientMatrix = fillSpecialInfluenceCoefficientMatrix(multipleCoef, elemLength, N, elementX)

    #2 write system.csv
    writeDataForSlau('system.csv', N, specialCoefficientMatrix, pressure)

    #3 exec slau
    execSlau('slau.exe')

    #4 read out.csv
    data.arrayD = readOut('out.csv')

    # count deformation energy
    deformationEnergy = 0.0
    for i in range(0, N):
        deformationEnergy = deformationEnergy - (data.arrayD[i] * pressure[i] * elemLength/2.0)
    return deformationEnergy



### math ################################################################
def countDerivativeDeformationEnergy(deformationEnergy1, deformationEnergy2, elemLength):
    return (deformationEnergy2 - deformationEnergy1) / (2 * elemLength)


def countIntensityCoefficient(derivativeDeformationEnergy, v, G):
    return math.sqrt(derivativeDeformationEnergy * G / (1 - v))


def countTheoreticalIntensityCoefficient(pressure, length,elemLength):
    return pressure[0] * math.sqrt(math.pi * (length - elemLength)/2 )
#    print "TheoreticalIntensityCoefficient = ", pressure[10] * math.sqrt(math.pi * (length - elemLength)/2 ) * (1 - 2.0/ math.pi * math.acos((0.2)))
#    return pressure[10] * math.sqrt(math.pi * (length - elemLength)/2 ) * (1 - 2.0/ math.pi * math.acos(0.2) )




def derivativeFxyy(X, V, Xi, Y, elemLength):
    return 1.0 / (4.0 * math.pi * (1 - V)) * (((X - Xi - elemLength / 2) ** 2 - Y ** 2) / (((X - Xi - elemLength / 2.0) ** 2 + Y ** 2) ** 2) - ((X - Xi + elemLength / 2) ** 2 - Y ** 2) / (((X - Xi + elemLength / 2) ** 2 + Y ** 2) ** 2))
    

def derivativeFyyy(X, V, Xi, Y, elemLength):
    return (2.0 * Y) / (4.0 * math.pi * (1 - V)) * ((X - Xi - elemLength / 2) / (((X - Xi - elemLength / 2) ** 2 + Y ** 2) ** 2) - (X - Xi + elemLength / 2) / (((X - Xi + elemLength / 2) ** 2 + Y ** 2) ** 2))


def derivativeFxx(X, V, Xi, Y, elemLength):
    return 1.0 / (4 * math.pi * (1 - V)) * ((X - Xi - elemLength / 2) / ((X - Xi - elemLength / 2) ** 2 + Y ** 2)  - (X - Xi + elemLength / 2) / ((X - Xi + elemLength / 2) ** 2 + Y ** 2))


def  derivativeFyy(X  , V  , Xi  , Y  , elemLength  ):
    return (-1.0) * derivativeFxx(X, V, Xi, Y, elemLength)


def  derivativeFxy(X  , V  , Xi  , Y  , elemLength  ):
    return 1.0 / (4 * math.pi * (1 - V)) * (Y / ((X - Xi - elemLength / 2) ** 2 + Y ** 2) - Y / ((X - Xi + elemLength / 2) ** 2 + Y ** 2))


def  derivativeFy(X  , V  , Xi  , Y  , elemLength  ):
    return (-1.0) / (4.0 * math.pi * (1 - V)) * (math.atan(Y / (X - Xi - elemLength / 2)) - math.atan(Y / (X - Xi + elemLength / 2)))


def  derivativeFx(X  , V  , Xi  , Y  , elemLength  ):
    return 1.0 / (4.0 * math.pi * (1 - V)) * (math.log((X - Xi - elemLength / 2) ** 2 + Y ** 2) - math.log((X - Xi + elemLength / 2) ** 2 + Y ** 2))


def  Uxi(Di  , X  , V  , Xi  , Y  , elemLength  ):
    return Di * ((-1.0) * (1 - 2.0 * V) * derivativeFx(X, V, Xi, Y, elemLength) - Y * derivativeFxy(X, V, Xi, Y, elemLength))


def  Uyi(Di  , X  , V  , Xi  , Y  , elemLength  ):
    return Di * (2.0 * (1 - V) * derivativeFy(X, V, Xi, Y, elemLength) - Y * derivativeFyy(X, V, Xi, Y, elemLength))


def  sigmaXXi(G  , Di  , X  , V  , Xi  , Y  , elemLength  ):
    return 2.0 * G * Di * (derivativeFyy(X, V, Xi, Y, elemLength) + Y * derivativeFyyy(X, V, Xi, Y, elemLength))


def  specialSigmaYYiForExternalOx(G , Di, X , V, elemLength, endOfFructure):
    r = abs (X - endOfFructure)
    return - G * Di / (2.0 * math.pi * (1 - V)) * 1.0 / math.sqrt(elemLength/2.0 * r) * ( math.atan(math.sqrt(elemLength/r )) - math.sqrt(elemLength * r) /(r + elemLength) )


def  sigmaYYi(G, Di  , X  , V  , Xi  , Y  , elemLength  ):
    return 2.0 * G * Di * (derivativeFyy(X, V, Xi, Y, elemLength) - Y * derivativeFyyy(X, V, Xi, Y, elemLength))


def  sigmaXYi(G, Di  , X  , V  , Xi  , Y  , elemLength  ):
    return 2.0 * G * Di * ((-1) * Y * derivativeFxyy(X, V, Xi, Y, elemLength))


def  sigmaXX(N   , G  , arrayD  , X  , V  , arrayX, Y  , elemLength  ):
    #print 'len(arrayD)',len(arrayD), 'len(arrayX)', len(arrayX)
    sigmaXX = 0.0
    for i in range(0, N):
        sigmaXX += sigmaXXi(G, arrayD[i], X, V, arrayX[i], Y, elemLength)
    return sigmaXX


def  sigmaYY(N   , G  , arrayD, X  , V  , arrayX, Y  , elemLength  ):
    sigmaYY = 0.0
    for i in range(0, N):
        sigmaYY += sigmaYYi(G, arrayD[i], X, V, arrayX[i], Y, elemLength)
    return sigmaYY


def  sigmaXY(N   , G  , arrayD, X  , V  , arrayX, Y  , elemLength  ):
    sigmaXY = 0.0
    for i in range(0, N):
        sigmaXY += sigmaXYi(G, arrayD[i], X, V, arrayX[i], Y, elemLength)
    return sigmaXY


def  specialSigmaYYforExternalOx( N, G  , arrayD, X  , V  , arrayX, elemLength, fromPoint, toPoint  ):
    theoreticalIntensityCoefficient = countTheoreticalIntensityCoefficient(data.pressure, data.length, data.elemLength)
    #print 'theoreticalIntensityCoefficient',theoreticalIntensityCoefficient
    sigmaTheoretical = countSigmaTheoretical (abs(X -fromPoint ), theoreticalIntensityCoefficient)

    assert ( not X>toPoint)
#    print "X=", X
    sigmaYY = 0.0
    for i in range(0, N):
        if (i != 0) and (i != N - 1):
            sigmaYY += sigmaYYi(G, arrayD[i], X, V, arrayX[i], 0, elemLength)
     #       print sigmaYY/sigmaTheoretical
        else:
            if i == 0:
                sigmaYY += specialSigmaYYiForExternalOx(G, arrayD[i], X, V, elemLength, fromPoint)
                #print "specialSigmaYYiForExternalOx i=0",specialSigmaYYiForExternalOx(G, arrayD[i], X, V, elemLength, fromPoint)/sigmaTheoretical, '\n'
     #           print sigmaYY/sigmaTheoretical
            if i == N - 1:
                sigmaYY += specialSigmaYYiForExternalOx(G, arrayD[i], X, V, elemLength, toPoint)
     #           print sigmaYY/sigmaTheoretical
    return sigmaYY


def countSigmaTheoretical (r, intensityCoefficient):
    return intensityCoefficient / math.sqrt(2 *math.pi * r)
#    return
#################################################################################




"""
def countArraySigmaYY(N, G, arrayD, v, elementX, elemLength, fromPoint, toPoint, offset, precisionCoefficient):
    arraySigmaYY = []
    X = 0.0
    Y = 0.0
    step = data.elemLength/precisionCoefficient
    assert(len(data.arrayD) != 0)
    for i in range(0, precisionCoefficient * data.N + 2* offset +1):
        arraySigmaYY.append([])
        for j in range(0, precisionCoefficient * data.N + 2 * offset + 1):
            X = data.fromPoint - offset * step + i  * step
            Y = data.fromPoint - offset * step + j  * step
            if((X > fromPoint - 1.2 * elemLength) and (X < toPoint + 1.2 * elemLength) and (Y > -1.2 * elemLength and Y < 1.2 * elemLength)):
                arraySigmaYY[i].insert(j, 0.0)
            else:
                arraySigmaYY[i].insert(j, sigmaYY(data.N, data.G, data.arrayD, X, data.v, data.elementX, Y, data.elemLength))
    return arraySigmaYY


def countArraySigmaXY(N, G, arrayD, v, elementX, elemLength, fromPoint, toPoint, offset, precisionCoefficient):
    arraySigmaXY = []
    X = 0.0
    Y = 0.0
    step = data.elemLength/precisionCoefficient
    assert(len(data.arrayD) != 0)
    for i in range(0, precisionCoefficient * data.N + 2* offset +1):
        arraySigmaXY.append([])
        for j in range(0, precisionCoefficient * data.N + 2 * offset + 1):
            X = data.fromPoint - offset * step + i  * step
            Y = data.fromPoint - offset * step + j  * step
            if((X > data.fromPoint - 1.2 * data.elemLength) and (X < data.toPoint + 1.2 * data.elemLength) and (Y > -1.2 * data.elemLength and Y < 1.2 * data.elemLength)):
                arraySigmaXY[i].insert(j, 0.0)
            else:
                arraySigmaXY[i].insert(j, sigmaXY(data.N, data.G, data.arrayD, X, data.v, data.elementX, Y, data.elemLength))
    return arraySigmaXY


def countArraySigmaXX(N, G, arrayD, v, elementX, elemLength, fromPoint, toPoint, offset, precisionCoefficient):
    arraySigmaXX = []
    X = 0.0
    Y = 0.0
    step = elemLength/precisionCoefficient
    assert(len(arrayD) != 0)
    for i in range(0, precisionCoefficient * N + 2* offset +1):
        arraySigmaXX.append([])
        for j in range(0, precisionCoefficient * N + 2 * offset + 1):
            X = fromPoint - offset * step + i  * step
            Y = fromPoint - offset * step + j  * step
            if((X > fromPoint - 1.2 * elemLength) and (X < toPoint + 1.2 * elemLength) and (Y > -1.2 * elemLength and Y < 1.2 * elemLength)):
                arraySigmaXX[i].insert(j, 0.0)
            else:
                arraySigmaXX[i].insert(j, sigmaXX(N, G, arrayD, X, v, elementX, Y, elemLength))
    return arraySigmaXX
"""


def countTableSigma(N, G, arrayD, v, elementX, elemLength, fromPoint, toPoint, offset, precisionCoefficient, charInfo):
    tableSigma =  TableSigma()
    tableSigma.charInfo = charInfo

    step = elemLength/precisionCoefficient
    assert(len(arrayD) != 0)
    xMax = precisionCoefficient * N + 2 * offset
    yMax = precisionCoefficient * N + 2 * offset

    for i in range(0, xMax +1):
        X = fromPoint - offset * step + i  * step
        tableSigma.X.insert(i, X)

    for j in range(0,  yMax+ 1):
        Y = fromPoint - offset * step + j  * step
        tableSigma.Y.insert(j, Y)

    print "xMax = ", xMax,"yMax = ", yMax
    print tableSigma.arraySigma

    additionlOffset = 1.2
    for i in range(0, xMax +1):
        tableSigma.arraySigma.append([])
        for j in range(0,  yMax+ 1):
            #print "i =", i
            #print "j =", j
            X = tableSigma.X[i]
            Y = tableSigma.Y[j]

            if(( X > fromPoint - additionlOffset * elemLength) and (X < toPoint + additionlOffset * elemLength) and (Y > - additionlOffset * elemLength and Y < additionlOffset * elemLength)):
                tableSigma.arraySigma[i].insert(j, 0.0)
            else:
                if tableSigma.charInfo == "sigmaXX" :
                    tableSigma.arraySigma[i].insert(j, sigmaXX(N, G, arrayD, X, v, elementX, Y, elemLength))

                if tableSigma.charInfo == "sigmaXY" :
                    tableSigma.arraySigma[i].insert(j, sigmaXY(N, G, arrayD, X, v, elementX, Y, elemLength))

                if tableSigma.charInfo == "sigmaYY" :
                    tableSigma.arraySigma[i].insert(j, sigmaYY(N, G, arrayD, X, v, elementX, Y, elemLength))

    return tableSigma



def plotSigma(tableSigma):
    import numpy as np
    import pylab as pl
    from matplotlib import cm
    import matplotlib
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    x = len(tableSigma.arraySigma)
    y = len(tableSigma.arraySigma[0])
    #X, Y = np.meshgrid(np.arange(-x/2.0, x/2.0, 1), np.arange(0,y,1))
    X, Y = np.meshgrid(tableSigma.Y, tableSigma.X)
    Z = np.zeros((x, y), 'Float32')
    for i in range(0, x):
        for j in range(0, y):
            Z[i,j] = tableSigma.arraySigma[i][j] / 10**6
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot_surface(X = X, Y = Y, Z = Z, color = 'GREEN', rstride=1, cstride=1, cmap=cm.jet , linewidth=0, antialiased = True)
    plt.title(tableSigma.charInfo + ", MPa")
    plt.xlabel('y, m')
    plt.ylabel('x, m')
    pl.savefig(tableSigma.charInfo)
    plt.show()
    return


def plotContour(tableSigma):
    arraySigma = tableSigma.arraySigma
    imageName = tableSigma.charInfo

    import numpy as np
    import pylab as pl
    from matplotlib import cm
    import matplotlib
    import matplotlib.pyplot as plt
    ## Option (True/False) to fill the gaps between contours :
    contourfillflag = True
    ## Option (True/False) to only fill the gaps between contours, without drawing
    ## the contours :
    contourfillonlyflag = True
    ## Contour levels can be True (automatic), the number of leves or the list of
    ## z levels :
    # contourlevels = True
    contourlevels = 100
    # contourlevels = linspace(0.0,1.0,10)
    # contourlevels = (0.0,0.2,0.4,0.6,0.8,1.0)
    ## Option (True/False) to show automatic contour labels  :
    contourlabelsflag = True
    ## Option (True/False) to show a vertical color bar of the contour levels  :
    contourlevelsbarflag = True
    ## Labels for x and y axis :
    xlabeltext = r'x'; ylabeltext = r'y'
    ## Plot title, here including TeX expressions (inside '$') :
    titletext = imageName
    ## Option (True/False) to show a grid of dashed lines :
    gridflag = True
    ##############################################################################
    x = len(arraySigma)
    y = len(arraySigma[0])
    X, Y = np.meshgrid(np.arange(-x/2.0,x/2.0,1), np.arange(0,y,1))
    Z = np.zeros((x,y), 'Float32')
    for i in range(0,x):
        for j in range(0,y):
            Z[i,j] = arraySigma[i][j]
    if contourfillflag:
        if contourlevels == True:
            if not contourfillonlyflag:
                c1 = pl.contour(x,y,z,colors='k')
            c2 = pl.contourf(X,Y,Z)
        else:
            if not contourfillonlyflag:
                c1 = pl.contour(X,Y,Z,contourlevels,colors='k')
            c2 = pl.contourf(X,Y,Z,contourlevels)
        if contourlevelsbarflag:
            cb2 = pl.colorbar(c2)
    else:
        if contourlevels == True:
            c1 = pl.contour(X,Y,Z)
        else:
            c1 = pl.contour(X,Y,Z,contourlevels)
        if contourlevelsbarflag:
            pl.colorbar(c1)
    if contourlabelsflag and ((not contourfillflag) or (not contourfillonlyflag)):
        l1 = clabel(c1)
    plt.xlabel(xlabeltext); plt.ylabel(ylabeltext); plt.title(titletext)
    plt.grid(gridflag)
    plt.savefig(imageName)
    plt.show()
    return



def countIntensityCoefficientByOffsetMethod(relationalOffset, N, G, arrayD, V, arrayX, elemLength, fromPoint, toPoint, length):
    X = fromPoint - relationalOffset * length/2.0
#    print "specialSigmaYYforExternal Ox = ", ( specialSigmaYYforExternalOx ( N, G, arrayD, X, V, arrayX, elemLength, fromPoint, toPoint) * math.sqrt(2 * math.pi * (fromPoint - X)) )
    return specialSigmaYYforExternalOx ( N, G, arrayD, X, V, arrayX, elemLength, fromPoint, toPoint) * math.sqrt(2 * math.pi * (fromPoint - X))


#function count table of K=K(x) for ploting
def countCorelationTableIntensityCoefficient(startDegree, finishDegree, numberOfPoints, pressure, N, G, arrayD, v, elementX, elemLength, fromPoint, toPoint, length):
    theoreticalIntensityCoefficient = countTheoreticalIntensityCoefficient(pressure, length, elemLength)
    print "theoreticalIntensityCoefficient =", theoreticalIntensityCoefficient

    table = CorelationTableIntensityCoefficient()
    stepDegree = (finishDegree - startDegree )*1.0 / numberOfPoints
    for i in range(0, numberOfPoints+1):
        table.relativeOffset.insert(i, (- (startDegree + i* stepDegree)))
        intensityCoefficientByOffsetMethod = countIntensityCoefficientByOffsetMethod(10 **(- (startDegree + i* stepDegree)), N, G, arrayD, v, elementX, elemLength,  fromPoint, toPoint, length)
        print "intensityCoefficientByOffsetMethod = ", intensityCoefficientByOffsetMethod
        table.relativeIntensityCoefficient.insert(i, intensityCoefficientByOffsetMethod / theoreticalIntensityCoefficient)

    return table


def plotTable(table, imageName):
    import matplotlib.pyplot as plt
    import pylab as pl
    plt.plot(table.relativeOffset, table.relativeIntensityCoefficient, linewidth = 3.0)
    plt.title(imageName)
    pl.savefig(imageName)
    plt.show()
    return



if __name__ == '__main__':
    #1 fill data
    data = readInput("input.csv")
    assert(int(data.N) == int(len(data.pressure)))
    data.length = abs(data.toPoint - data.fromPoint)
    data.elemLength = data.length / data.N
    data.multipleCoef = (-data.G) / (math.pi * (1 - data.v))
    data.elementX = fillElementX(data.N, data.elemLength, data.fromPoint)

    #2 count deformation energy with initial data
    deformationEnergy1 = countDeformationEnergy(data.multipleCoef, data.elemLength, data.N, data.elementX, data.pressure)


    #2.1 count sigmaXX and sigmaYY  and   plot sigmaXX and sigmaYY
    charInfo = "sigmaXX"
    tableSigma = countTableSigma(data.N, data.G, data.arrayD, data.v, data.elementX, data.elemLength, data.fromPoint, data.toPoint, 20, 1, charInfo)
    plotSigma(tableSigma)
    plotContour(tableSigma)


    charInfo = "sigmaYY"
    tableSigma = countTableSigma(data.N, data.G, data.arrayD, data.v, data.elementX, data.elemLength, data.fromPoint, data.toPoint, 20, 1, charInfo)
    plotSigma(tableSigma)
    plotContour(tableSigma)


    charInfo = "sigmaXY"
    tableSigma = countTableSigma(data.N, data.G, data.arrayD, data.v, data.elementX, data.elemLength, data.fromPoint, data.toPoint, 20, 1, charInfo)
    plotSigma(tableSigma)
    plotContour(tableSigma)

    #import matplotlib.pyplot as plt
    #plt.show()

    #2.2.1
    #theoreticalIntensityCoefficient = countTheoreticalIntensityCoefficient(data.pressure, data.length, data.elemLength)
    #print 'theoreticalIntensityCoefficient',theoreticalIntensityCoefficient

    #2.3 plotting the relative stress intensity factor of the relative    #count corelation intensity ciefficient
    table = countCorelationTableIntensityCoefficient(1, 7, 10, data.pressure, data.N, data.G, data.arrayD, data.v, data.elementX, data.elemLength,  data.fromPoint, data.toPoint, data.length)


    #2.4 plot relativeOffset and relativeIntensityCoefficient
    plotTable(table, 'relativeIntensityCoefficient')

    """
    import numpy as np
    import pylab as pl
    #from numpy import *
    from matplotlib import cm
    import matplotlib
    import matplotlib.pyplot as plt2
    from mpl_toolkits.mplot3d import Axes3D


    #fig = plt.figure()
    #ax = Axes3D(fig)



    plt2.plot(table.relativeOffset, table.relativeIntensityCoefficient)
    plt2.show()

    #plt.imsave ("qwer")
    #plt.savefig("plot")

    #plt.show()




    #compare with theoretical intensity coefficient
    theoreticalIntensityCoefficient = countTheoreticalIntensityCoefficient(data.pressure, data.length, data.elemLength)
    intensityCoefficientByOffsetMethod = countIntensityCoefficientByOffsetMethod(0.01, data.N, data.G, data.arrayD, data.v, data.elementX, data.elemLength,  data.fromPoint, data.toPoint, data.length)
    # print 'intensityCoefficientByOffsetMethod1 = ', intensityCoefficientByOffsetMethod, '\t  ', 0.01,'\t  ', intensityCoefficientByOffsetMethod/ theoreticalIntensityCoefficient

    intensityCoefficientByOffsetMethod = countIntensityCoefficientByOffsetMethod(0.0001, data.N, data.G, data.arrayD, data.v, data.elementX, data.elemLength, data.fromPoint, data.toPoint,  data.length)
    # print 'intensityCoefficientByOffsetMethod2 = ', intensityCoefficientByOffsetMethod, '\t  ', 0.0001,'\t  ', intensityCoefficientByOffsetMethod/ theoreticalIntensityCoefficient

    intensityCoefficientByOffsetMethod = countIntensityCoefficientByOffsetMethod(0.000001, data.N, data.G, data.arrayD, data.v, data.elementX, data.elemLength, data.fromPoint, data.toPoint, data.length)
    # print 'intensityCoefficientByOffsetMethod3 = ', intensityCoefficientByOffsetMethod, '\t  ', 0.000001,'\t  ', intensityCoefficientByOffsetMethod/ theoreticalIntensityCoefficient


    """

    """
    #3 change data
    data.fromPoint -= data.elemLength
    data.toPoint += data.elemLength
    data.N += 2
    data.length = abs(data.toPoint - data.fromPoint)
    data.elementX = fillElementX(data.N, data.elemLength, data.fromPoint)
    for i in range(1, len(data.pressure)):
        assert(data.pressure[i] == data.pressure[i - 1])

    data.pressure.insert(0, data.pressure[0])
    data.pressure.insert(len(data.pressure), data.pressure[0])


    #4 count deformation energy with changed data
    deformationEnergy2 = countDeformationEnergy(data.multipleCoef, data.elemLength, data.N, data.elementX, data.pressure)


    #5 count deformation energy derivative
    derivativeDeformationEnergy = countDerivativeDeformationEnergy(deformationEnergy1, deformationEnergy2, data.elemLength)
    print 'derivativeDeformationEnergy',derivativeDeformationEnergy

    #6 count intensity coefficient
    intensityCoefficient = countIntensityCoefficient(derivativeDeformationEnergy, data.v, data.G)
    print 'intensityCoefficient',intensityCoefficient

    theoreticalIntensityCoefficient = countTheoreticalIntensityCoefficient(data.pressure, data.length, data.elemLength)
    print 'theoreticalIntensityCoefficient',theoreticalIntensityCoefficient
    print '\n', intensityCoefficient/theoreticalIntensityCoefficient
    """