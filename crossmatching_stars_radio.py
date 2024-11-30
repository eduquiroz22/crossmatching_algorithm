import numpy as np
import pandas as pd
from datetime import datetime
from multiprocessing import Pool
from functools import partial

directory = '../0 survey/'
VLASS_file =  directory+'first.csv'	
EDR3_file = directory+'4mas_dr3_predicted2.csv'					
name_output = 'matches_first_2sec_250pc_3.csv'
max_dist = 2 #arcsec
max_dif = max_dist/3600
n_proc = 6

#Gaia: RA1, Dec1
#VLASS: RA2, Dec2

def calculateDistance(RA_1,RA_2,Dec_1,Dec_2):
    dRA = abs(RA_2-RA_1)
    y = np.sqrt((np.cos(np.radians(Dec_2))*np.sin(np.radians(dRA)))**2+(np.cos(np.radians(Dec_1))*np.sin(np.radians(Dec_2))-np.sin(np.radians(Dec_1))*np.cos(np.radians(Dec_2))*np.cos(np.radians(dRA)))**2)
    x = np.sin(np.radians(Dec_1))*np.sin(np.radians(Dec_2))+np.cos(np.radians(Dec_1))*np.cos(np.radians(Dec_2))*np.cos(np.radians(dRA))
    dist = np.arctan2(y,x)
    dist = np.degrees(dist)*3600.
    return dist

def diference_ra(ra1, ra2):
    dif = abs(ra2-ra1)
    if dif>180:
        dif = 360-dif
    return dif


def binary_search_ra(arr, x, dif=max_dif, mid=0, low=0, high=0):
    if high==0:
        high=len(arr)-1
    result = []
 
    while low <= high:
        mid = round((high + low) / 2)
 
        # If x is greater, ignore left half
        if arr[mid]-x < -dif:
            low = mid + 1
 
        # If x is smaller, ignore right half
        elif arr[mid]-x > dif:
            high = mid - 1
 
        # means x is present at mid
        else:
            c=1
            while diference_ra(arr[mid-c],x)<=dif:
                c+=1
            result.append(mid-c+1)

            l = len(arr)
            c = 1
            while diference_ra(arr[(mid+c)%l],x)<=dif:
                c+=1
            result.append((mid+c-1)%l)
            return result
 
    # If we reach here, then the element was not present
    return result


def binary_search(arr, x, dif=max_dif, mid=0, low=0, high=0):
    if high==0:
        high=len(arr)-1
    result = []
 
    while low <= high:
        mid = round((high + low) / 2)
 
        # If x is greater, ignore left half
        if arr[mid]-x < -dif:
            low = mid + 1
 
        # If x is smaller, ignore right half
        elif arr[mid]-x > dif:
            high = mid - 1
 
        # means x is present at mid
        else:
            if mid>0:
                c=1
                while abs(arr[mid-c]-x)<=dif:
                    c+=1
                    if (mid-c)<0:
                        break
                result.append(mid-c+1)
            else:
                result.append(mid)
            if mid<len(arr)-1:
                c = 1
                while abs(arr[mid+c]-x)<=dif:
                    c+=1
                    if (mid+c)>=len(arr):
                        break
                result.append(mid+c-1)
            else:
                result.append(mid)
            return result
 
    # If we reach here, then the element was not present
    return result


#(1,5,7,9)
#(0,1,2,3)
start = datetime.now() 
print("leyendo archivo gaia...")
source_id, RA1_, Dec1_, parallax, pmra, pmdec, epoch1, RA1, Dec1 = np.loadtxt(EDR3_file,delimiter=',',usecols=(0,1,2,3,4,5,6,7,8),unpack=True, skiprows=1)
source_id = np.loadtxt(EDR3_file,delimiter=',',usecols=(0),unpack=True, skiprows=1, dtype=str)
end = datetime.now()
print(f'Tiempo gaia: {end-start}')

start = datetime.now() 
print("leyendo archivo first...")
RA2, Dec2, epoch2 = np.loadtxt(VLASS_file,delimiter=',',usecols=(0,1,2),unpack=True, skiprows=1)
end = datetime.now()
print(f'Tiempo vlass: {end-start}')

print("preparando data...")

index1 = np.arange(len(RA1))
index2 = np.arange(len(RA2))

gaia = [index1, Dec1, RA1, parallax, source_id, epoch1, Dec1_, RA1_, pmra, pmdec] 
vlass = [index2, Dec2, RA2, epoch2]

order_Dec1 = np.argsort(Dec1)
order_Dec2 = np.argsort(Dec2)

gaia = [e[order_Dec1] for e in gaia]  
vlass = [e[order_Dec2] for e in vlass]

min_dec = (gaia[1][0],1) if gaia[1][0]>vlass[1][0] else (vlass[1][0],2)
max_dec = (gaia[1][-1],1) if gaia[1][-1]<vlass[1][-1] else (vlass[1][-1],2)

mi = 0
while True:
    min_dec = (gaia[1][0+mi],1) if gaia[1][0]>vlass[1][0] else (vlass[1][0+mi],2)
    try:
        arr = vlass[1] if min_dec[1]==1 else gaia[1]
        min_index = binary_search(arr, min_dec[0])[0]
        break
    except:
        mi+=1

ma = 0
while True:
    max_dec = (gaia[1][-1-ma],1) if gaia[1][-1]<vlass[1][-1] else (vlass[1][-1-ma],2)
    try:
        arr = vlass[1] if max_dec[1]==1 else gaia[1]
        max_index = binary_search(arr, max_dec[0])[1]
        break
    except:
        ma+=1


if min_dec[1]==1:
    if max_dec[1]==1:
        vlass = [e[min_index:max_index+1] for e in vlass]
    else:
        vlass = [e[min_index:] for e in vlass]
        gaia = [e[:max_index+1] for e in gaia]
else:
    if max_dec[1]==1:
        vlass = [e[:max_index+1] for e in vlass]
        gaia = [e[min_index:] for e in gaia]
    else:
        gaia = [e[min_index:max_index+1] for e in gaia]


def get_max_dif_ra(x, max_dif):
    x = np.radians(x)
    return np.degrees(np.arccos((np.cos(x))**(-2)*np.cos(np.radians(max_dif))-(np.tan(x))**2))

def ra_iteration(result, l):
    ini = result[0]
    fin = result[1]

    if ini<=fin:
        iteration = [i for i in range(ini,fin+1)]
    else:
        arr_it = [i for i in range(ini,l)] + [i for i in range(0,ini+1)]
    return iteration


def comparar(n):
    # Test array
    matches = []
    arr = vlass[1]
    dec = gaia[1][n]
    
    # Function call
    resultDec = binary_search(arr, dec)
    
    if len(resultDec)!=0:
        # Test array
        vlass_copy = [e[resultDec[0]:resultDec[-1]+1].copy() for e in vlass]
        order = np.argsort(vlass_copy[2])
        vlass_copy = [e[order] for e in vlass_copy]
        arr = vlass_copy[2]
        ra = gaia[2][n]

        max_dif_ra = get_max_dif_ra(abs(dec)+max_dif,max_dif)
        # Function call
        if (ra<max_dif_ra or ra>(360-max_dif_ra)) or ((arr[0]-max_dif_ra)<ra and ra<(arr[-1]+max_dif_ra)):
            result = binary_search_ra(arr, ra, max_dif_ra)
            r = []
            if len(result)!=0:
                for i in ra_iteration(result, len(arr)):
                    dist = calculateDistance(ra, arr[i], dec, vlass_copy[1][i])
                    if dist<max_dist and vlass_copy[3][i]!=2.0:
                        r.append([dist, i])
                
                for e in r:
                    i = e[1]
                    dist = e[0]
                    #print(n, "gaia: ", gaia[0][n], " vlass: ", vlass[0][i])
                    matches.append([gaia[4][n], gaia[3][n], gaia[0][n], gaia[7][n], gaia[6][n], vlass_copy[0][i], vlass_copy[2][i], vlass_copy[1][i], gaia[5][n], vlass_copy[3][i], dist, gaia[8][n], gaia[9][n]])
                                #['source_id', 'parallax', 'gaia_index','gaia_ra','gaia_dec', 'first_index', 'first_ra', 'first_dec', 'gaia_epoch', 'first_epoch', 'distance', 'pmra', 'pmdec'])
#gaia = [index1, Dec1, RA1, parallax, source_id, epoch1, Dec1_, RA1_, pmra, pmdec] 
#vlass = [index2, Dec2, RA2, epoch2]
                return matches
    return 0

print("empezando comparacion...")


if __name__  == '__main__':
    longitude = len(gaia[0])
    index_g = [i for i in range(longitude)]
    start = datetime.now()
    with Pool(processes=n_proc) as pool:
        matches_ = pool.map(comparar, index_g)
    
    matches = []
    for e in matches_:
    	if e!=0:
    		matches = matches + e
    
    print("Coincidencias hayadas: ", len(matches))    
    end = datetime.now()
    
    print("terminado")
    df = pd.DataFrame(matches, columns = ['source_id', 'parallax', 'gaia_index','gaia_ra','gaia_dec', 'first_index', 'first_ra', 'first_dec', 'gaia_epoch', 'first_epoch', 'distance', 'pmra', 'pmdec'])
    df = df.sort_values(by='parallax').reset_index(drop=True)
    df.to_csv(name_output,index=False)
    print(f'Tiempo: {end-start}')
