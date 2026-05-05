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
    data=False
    nlines1=0
    for line in file:
        if line.strip()=="*BEGIN_DATA":
            data=True
            continue
        elif line.strip()=="*END":
            data=False
            continue
        elif data:
            nlines1 += 1
    wavelength1 = np.zeros(nlines1)
    n1 = np.zeros(nlines1)
    k1 = np.zeros(nlines1)

with open(filename, 'r') as file:
    data=False
    i = 0
    for line in file:
        if line.strip()=="*BEGIN_DATA":
            data=True
            continue
        elif line.strip()=="*END":
            data=False
            continue
        elif data:
            line = line.strip()
            columns = line.split()
            wavelength1[i] = float(columns[0])
            n1[i] = float(columns[1])
            k1[i] = float(columns[2])
            i += 1

if (len(sys.argv) > 2):
    with open(filename2, 'r') as file:
        data=False
        nlines2=0
        for line in file:
            if line.strip()=="*BEGIN_DATA":
                data=True
                continue
            elif line.strip()=="*END":
                data=False
                continue
            elif data:
                nlines2 += 1
        wavelength2 = np.zeros(nlines2)
        n2 = np.zeros(nlines2)
        k2 = np.zeros(nlines2)
    
    with open(filename2, 'r') as file:
        data=False
        i = 0
        for line in file:
            if line.strip()=="*BEGIN_DATA":
                data=True
                continue
            elif line.strip()=="*END":
                data=False
                continue
            elif data:
                line = line.strip()
                columns = line.split()
                wavelength2[i] = float(columns[0])
                n2[i] = float(columns[1])
                k2[i] = float(columns[2])
                i += 1

if (len(sys.argv) > 3):
    with open(filename3, 'r') as file:
        data=False
        nlines3=0
        for line in file:
            if line.strip()=="*BEGIN_DATA":
                data=True
                continue
            elif line.strip()=="*END":
                data=False
                continue
            elif data:
                nlines3 += 1
        wavelength3 = np.zeros(nlines3)
        n3 = np.zeros(nlines3)
        k3 = np.zeros(nlines3)
    
    with open(filename3, 'r') as file:
        data=False
        i = 0
        for line in file:
            if line.strip()=="*BEGIN_DATA":
                data=True
                continue
            elif line.strip()=="*END":
                data=False
                continue
            elif data:
                line = line.strip()
                columns = line.split()
                wavelength3[i] = float(columns[0])
                n3[i] = float(columns[1])
                k3[i] = float(columns[2])
                i += 1

if (len(sys.argv) > 4):
    with open(filename4, 'r') as file:
        data=False
        nlines4=0
        for line in file:
            if line.strip()=="*BEGIN_DATA":
                data=True
                continue
            elif line.strip()=="*END":
                data=False
                continue
            elif data:
                nlines4 += 1
        wavelength4 = np.zeros(nlines4)
        n4 = np.zeros(nlines4)
        k4 = np.zeros(nlines4)
    
    with open(filename4, 'r') as file:
        data=False
        i = 0
        for line in file:
            if line.strip()=="*BEGIN_DATA":
                data=True
                continue
            elif line.strip()=="*END":
                data=False
                continue
            elif data:
                line = line.strip()
                columns = line.split()
                wavelength4[i] = float(columns[0])
                n4[i] = float(columns[1])
                k4[i] = float(columns[2])
                i += 1

if (len(sys.argv) > 5):
    with open(filename5, 'r') as file:
        data=False
        nlines5=0
        for line in file:
            if line.strip()=="*BEGIN_DATA":
                data=True
                continue
            elif line.strip()=="*END":
                data=False
                continue
            elif data:
                nlines5 += 1
        wavelength5 = np.zeros(nlines5)
        n5 = np.zeros(nlines5)
        k5 = np.zeros(nlines5)
    
    with open(filename5, 'r') as file:
        data=False
        i = 0
        for line in file:
            if line.strip()=="*BEGIN_DATA":
                data=True
                continue
            elif line.strip()=="*END":
                data=False
                continue
            elif data:
                line = line.strip()
                columns = line.split()
                wavelength5[i] = float(columns[0])
                n5[i] = float(columns[1])
                k5[i] = float(columns[2])
                i += 1

if (len(sys.argv) > 6):
    with open(filename6, 'r') as file:
        data=False
        nlines6=0
        for line in file:
            if line.strip()=="*BEGIN_DATA":
                data=True
                continue
            elif line.strip()=="*END":
                data=False
                continue
            elif data:
                nlines6 += 1
        wavelength6 = np.zeros(nlines6)
        n6 = np.zeros(nlines6)
        k6 = np.zeros(nlines6)
    
    with open(filename6, 'r') as file:
        data=False
        i = 0
        for line in file:
            if line.strip()=="*BEGIN_DATA":
                data=True
                continue
            elif line.strip()=="*END":
                data=False
                continue
            elif data:
                line = line.strip()
                columns = line.split()
                wavelength6[i] = float(columns[0])
                n6[i] = float(columns[1])
                k6[i] = float(columns[2])
                i += 1

fig=plt.figure()
ax1 = fig.add_subplot(121)
ax1.plot(wavelength1,n1, color='blue', label=filename)
if (len(sys.argv) > 2):
    ax1.plot(wavelength2,n2, color='green', label=filename2)
if (len(sys.argv) > 3):
    ax1.plot(wavelength3,n3, color='red', label=filename3)
if (len(sys.argv) > 4):
    ax1.plot(wavelength4,n4, color='cyan', label=filename4)
if (len(sys.argv) > 5):
    ax1.plot(wavelength5,n5, color='magenta', label=filename5)
if (len(sys.argv) > 6):
    ax1.plot(wavelength6,n6, color='orange', label=filename6)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_title('n')
plt.legend()
ax2 = fig.add_subplot(122)
ax2.plot(wavelength1,k1, color='blue')
if (len(sys.argv) > 2):
    ax2.plot(wavelength2,k2, color='green')
if (len(sys.argv) > 3):
    ax2.plot(wavelength3,k3, color='red')
if (len(sys.argv) > 4):
    ax2.plot(wavelength4,k4, color='cyan')
if (len(sys.argv) > 5):
    ax2.plot(wavelength5,k5, color='magenta')
if (len(sys.argv) > 6):
    ax2.plot(wavelength6,k6, color='orange')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_title('k')
plt.tight_layout()
plt.show()
