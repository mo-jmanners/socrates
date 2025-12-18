import numpy as np
import sys, math
import matplotlib.pyplot as plt

if __name__ == '__main__':
    if (len(sys.argv) > 6):
        filename = sys.argv[1]
        filename2 = sys.argv[2]
        filename3 = sys.argv[3]
        filename4 = sys.argv[4]
        filename5 = sys.argv[5]
        filename6 = sys.argv[6]
    elif (len(sys.argv) > 5):
        filename = sys.argv[1]
        filename2 = sys.argv[2]
        filename3 = sys.argv[3]
        filename4 = sys.argv[4]
        filename5 = sys.argv[5]
    elif (len(sys.argv) > 4):
        filename = sys.argv[1]
        filename2 = sys.argv[2]
        filename3 = sys.argv[3]
        filename4 = sys.argv[4]
    elif (len(sys.argv) > 3):
        filename = sys.argv[1]
        filename2 = sys.argv[2]
        filename3 = sys.argv[3]
    elif (len(sys.argv) > 2):
        filename = sys.argv[1]
        filename2 = sys.argv[2]
    elif (len(sys.argv) > 1):
        filename = sys.argv[1]
    else:
        raise RuntimeError('Please enter .dat files to plot')

with open(filename, 'r') as file:
    nlines1 = sum(1 for _ in file)
    wavelength1 = np.zeros(nlines1)
    xsc1 = np.zeros(nlines1)

with open(filename, 'r') as file:
    i = 0
    for line in file:
        line = line.strip()
        columns = line.split()
        wavelength1[i] = float(columns[0])
        xsc1[i] = float(columns[1])
        i += 1

if (len(sys.argv) > 2):
    with open(filename2, 'r') as file:
        nlines2 = sum(1 for _ in file)
        wavelength2 = np.zeros(nlines2)
        xsc2 = np.zeros(nlines2)

    with open(filename2, 'r') as file:
        i = 0
        for line in file:
            line = line.strip()
            columns = line.split()
            wavelength2[i] = float(columns[0])
            xsc2[i] = float(columns[1])
            i += 1

if (len(sys.argv) > 3):
    with open(filename3, 'r') as file:
        nlines3 = sum(1 for _ in file)
        wavelength3 = np.zeros(nlines3)
        xsc3 = np.zeros(nlines3)

    with open(filename3, 'r') as file:
        i = 0
        for line in file:
            line = line.strip()
            columns = line.split()
            wavelength3[i] = float(columns[0])
            xsc3[i] = float(columns[1])
            i += 1

if (len(sys.argv) > 4):
    with open(filename4, 'r') as file:
        nlines4 = sum(1 for _ in file)
        wavelength4 = np.zeros(nlines4)
        xsc4 = np.zeros(nlines4)

    with open(filename4, 'r') as file:
        i = 0
        for line in file:
            line = line.strip()
            columns = line.split()
            wavelength4[i] = float(columns[0])
            xsc4[i] = float(columns[1])
            i += 1

if (len(sys.argv) > 5):
    with open(filename5, 'r') as file:
        nlines5 = sum(1 for _ in file)
        wavelength5 = np.zeros(nlines5)
        xsc5 = np.zeros(nlines5)

    with open(filename5, 'r') as file:
        i = 0
        for line in file:
            line = line.strip()
            columns = line.split()
            wavelength5[i] = float(columns[0])
            xsc5[i] = float(columns[1])
            i += 1

if (len(sys.argv) > 6):
    with open(filename6, 'r') as file:
        nlines6 = sum(1 for _ in file)
        wavelength6 = np.zeros(nlines6)
        xsc6 = np.zeros(nlines6)

    with open(filename6, 'r') as file:
        i = 0
        for line in file:
            line = line.strip()
            columns = line.split()
            wavelength6[i] = float(columns[0])
            xsc6[i] = float(columns[1])
            i += 1

fig=plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(wavelength1,xsc1, color='blue')
if (len(sys.argv) > 2):
    ax1.plot(wavelength2,xsc2, color='green')
if (len(sys.argv) > 3):
    ax1.plot(wavelength3,xsc3, color='red')
if (len(sys.argv) > 4):
    ax1.plot(wavelength4,xsc4, color='cyan')
if (len(sys.argv) > 5):
    ax1.plot(wavelength5,xsc5, color='magenta')
if (len(sys.argv) > 6):
    ax1.plot(wavelength6,xsc6, color='orange')
ax1.set_yscale('log')
plt.tight_layout()
plt.show()
