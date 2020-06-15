import pandas as pd
import numpy as np
import math
import os
from os import listdir
from os.path import isfile, join
from sklearn.metrics.pairwise import euclidean_distances
import time
import gzip



def getListOfFiles(dirName):
    # create a list of file and sub directories
    # names in the given directory
    listOfFile = os.listdir(dirName)
    allFiles = list()
    # Iterate over all the entries
    for entry in listOfFile:
        # Create full path
        fullPath = os.path.join(dirName, entry)
        # If entry is a directory then get the list of files in this directory
        if os.path.isdir(fullPath):
            allFiles = allFiles + getListOfFiles(fullPath)
        else:
            allFiles.append(fullPath)

    return allFiles

path=os.getcwd()+'/base'
pathoutput=os.getcwd()+'/output'
print(path,pathoutput)
outputfiles=getListOfFiles(pathoutput)
files=getListOfFiles(path)
outputfiles=[os.path.split(x)[1]  for x in outputfiles]
AROMTRP= ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2']
AROMPHE= ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']
AROMTYR= ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']
for pdbfile in files:
    print(pdbfile)
    if(os.path.split(pdbfile)[1]+'.txt' in outputfiles):
        print('Results already on database')
        continue
    #PDB = open(pdbfile, 'r')
    #AA=['GLY', 'ALA', 'LEU', 'VAL', 'ILE', 'PRO', 'PHE', 'SER', 'THR', 'CYS', 'TYR', 'ASN', 'GLN', 'ASP', 'BASP', 'AASP,' 'GLU', 'ARG', 'LYS', 'HIS', 'TRP', 'MET', 'MG', 'CU', 'FE2', 'FE', 'NI', 'CL', 'NA', 'MO1', 'MO3', 'MO4','MO5','MO6','MO7','MO8','MO9', 'CA', ' I', 'BR', 'K']

    datapdb=pd.DataFrame(columns=['Type','Atom ID','Atom Type','aa','Chain ID','X','Y','Z','occupancy','temperature factor','element symbol'])
    i=0
    print("Starting to load data")
    start_time = time.time()
    Amin=''
    Aromaticpos=[]
    AromaticArray={}
    AromaticPoints=[]
    AromaticNormals={}
    Invalids=[]
    if(pdbfile[-3:]=='.gz'):
        f=gzip.open(pdbfile, 'rt')
    elif(pdbfile[-4:]=='.pdb'):
        f=open(pdbfile, 'r')
    for line in f:
        if("ENDMDL" in line):
            break
        if (line[16:20].strip() in "HOH"):
            continue
        if("ATOM" in line [:6] or "HETATM" in line [:6]):

            if (line[17:20].strip() in ['TYR','PHE','TRP']):
                #print("{}{}".format(Amin,line[17:20].strip()))
                if(Amin==''):
                    Amin=line[20:27].strip()
                elif(Amin!=line[20:27].strip()):
                    AromaticArray[Amin]=np.asarray(Aromaticpos)
                    if(len(AromaticPoints)==3):
                        veca=np.subtract(AromaticPoints[1], AromaticPoints[0])
                        vecb=np.subtract(AromaticPoints[2], AromaticPoints[0])
                        AromaticNormals[Amin]=np.cross(veca, vecb)
                    else:
                        print("Incomplete Aromatic Structure at {}".format(Amin))
                        Invalids.append(Amin)
                    Amin=line[20:27].strip()
                    Aromaticpos=[]
                    AromaticPoints=[]
                elif((line[17:20].strip()=='TYR' and line[11:16].strip() in AROMTYR) or (line[17:20].strip()=='PHE' and line[11:16].strip() in AROMPHE) or (line[17:20].strip()=='TRP' and line[11:16].strip() in AROMTRP)):
                    if(Aromaticpos==[]):
                        Aromaticpos=[float(line[27:38].strip()),float(line[38:46].strip()),float(line[46:54].strip())]
                    else:
                        Aromaticpos=[(x+y)/2 for x, y in zip(Aromaticpos, [float(line[27:38].strip()),float(line[38:46].strip()),float(line[46:54].strip())])]
                    if(len(AromaticPoints)<3):
                        AromaticPoints.append([float(line[27:38].strip()),float(line[38:46].strip()),float(line[46:54].strip())])
            datapdb.loc[i] = [line[:6].strip(),line[6:11].strip(),line[11:16].strip(),line[17:20].strip(),line[20:27].strip(),line[27:38].strip(),line[38:46].strip(),line[46:54].strip(),line[54:60].strip(),line[60:66].strip(),line[66:80].strip()]
            #datapdb.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
            i+=1
    AromaticArray[Amin]=np.asarray(Aromaticpos)
    if(len(AromaticPoints)==3):
        veca=np.subtract(AromaticPoints[1], AromaticPoints[0])
        vecb=np.subtract(AromaticPoints[2], AromaticPoints[0])
        AromaticNormals[Amin]=np.cross(veca, vecb)
    elif(0<len(AromaticPoints)<3):
        print("Incomplete Aromatic Structure at {}".format(Amin))
        Invalids.append(Amin)
    print("DataLoaded")
    print("---%s seconds ---" % (time.time() - start_time))
    if 'Atom ID' in datapdb.columns:
        datapdb=datapdb.drop(['Atom ID','occupancy','temperature factor','element symbol'], axis=1)
    datapdb["X"] = np.float32(datapdb["X"], downcast='signed')
    datapdb["Y"] = np.float32(datapdb["Y"], downcast='signed')
    datapdb["Z"] = np.float32(datapdb["Z"], downcast='signed')

    dist = euclidean_distances(np.float32(datapdb[["X","Y","Z"]].to_numpy()), np.float32(datapdb[["X","Y","Z"]].to_numpy()))
    Dist= pd.DataFrame(data=dist, index=None)

    New=pd.concat([datapdb, Dist], axis=1, sort=False)

    hb=0
    lighbacep=['OD1', 'OD2', 'OE1', 'OE2', 'OG', 'OG1'] #3
    lighbdono= ['HE1', 'H', 'HE', 'HH11', 'HH12', 'HH21', 'HH22', 'HD1', 'HE2', 'HG', 'HG1', 'HG21', 'HG22', 'HG23', 'HH', 'HD21', 'HD22', 'HE21', 'HE22'] #18

    sb=0
    ligsb1=['OD2', 'OE2', 'OH']
    ligsb2=['NZ', 'NH2', 'NH1']
    aasbneg=['ASP', 'GLU']
    aasbpos= ['LYS', 'ARG']

    db=0
    ligdb=['SG']

    vdw=0
    ligvdw=['CB', 'CG1', 'CG2', 'CD1', 'CD2', 'CE']
    aavdw=['VAL', 'TRE', 'MET', 'LEU', 'ILE']

    lpi=0
    tshaped=0
    inter=0
    paralel=0
    aapi=['TYR','PHE', 'TRP']

    ctn=0
    ligctn=['MG', 'CU', 'K', 'FE2', 'FE', 'NI', 'NA', 'MO1', 'MO3', 'MO4','MO5','MO6','MO7','MO8','MO9', 'NZ', 'NH2', 'NH1']
    ligctn2=['CG','CE2','CG']
    aactn=['TYR','PHE', 'TRP']

    an=0
    ligan1=['CL', 'BR', 'I','OD2', 'OE2', 'OH']
    ligan2=['CG','CE2','CG']
    aaan=['TYR','PHE', 'TRP']

    spi=0
    ligspi1=['SG']
    ligspi2=['CG','CE2','CG']
    aaspi=['TYR','PHE', 'TRP']

    f= open('output/'+os.path.split(pdbfile)[1] +"-ext"+".txt","w+")
    print("Starting to process data")
    start_time = time.time()
    Exclusions=[]
    for i in range(len(New)):
        for j in range(i+1, len(New)):
            distance = New[j].iloc[i]
            if(distance>8 or distance ==0):
                continue
            atom1=New['Atom Type'].iloc[i]
            atom2=New['Atom Type'].iloc[j]
            aa1=New['aa'].iloc[i]
            aa2=New['aa'].iloc[j]
            chaincode1=New['Chain ID'].iloc[i]
            chaincode2=New['Chain ID'].iloc[j]
            distance=New[j].iloc[i]
            #print(chaincode1,chaincode2)
            #HYDROGEN BOND---HYDROGEN BOND---HYDROGEN BOND---HYDROGEN BOND---HYDROGEN BOND---HYDROGEN BOND---
            if atom1 in lighbacep[:] and atom2 in lighbdono[:] and 0.0 < distance < 3.1:
                hb+=1
                f.write('Hydrogen_Bond' +'\t\t'+ atom1 +'\t\t'+ aa1 +'\t\t'+ chaincode1 +'\t\t'+ atom2 +'\t\t'+ aa2 +'\t\t'+ chaincode2 +'\t\t' + str(distance)+ '\n')
            elif atom1 in lighbdono[:] and atom2 in lighbacep[:] and 0.0 < distance < 3.1:
                hb+=1
                f.write('Hydrogen_Bond' + '\t\t' + atom1 + '\t\t' + aa1 + '\t\t' + chaincode1 + '\t\t'+ atom2 + '\t\t' + aa2 + '\t\t' + chaincode2 +'\t\t' + str(distance)+ '\n')
            #SALTY_BRIDGE-----SALTY_BRIDGE-----SALTY_BRIDGE-----SALTY_BRIDGE-----SALTY_BRIDGE-----SALTY_BRIDGE-----
            if aa1 in aasbpos[:] and atom1 in ligsb2[:] and aa2 in aasbneg[:] and atom2 in ligsb1[:] and 0.0 < distance < 4.0:
                sb+=1
                f.write('Salt_Bridge'+ '\t\t' + atom1+ '\t\t' + aa1+ '\t\t' + chaincode1+ '\t\t' + atom2+ '\t\t' + aa2+ '\t\t' + chaincode2+ '\t\t' + str(distance)+ '\n')
            elif aa1 in aasbneg[:] and atom1 in ligsb1[:] and aa2 in aasbpos[:] and atom2 in ligsb2[:] and 0.0 < distance < 5.0:
                sb+=1
                f.write('Salt_Bridge' + '\t\t' + atom1 + '\t\t' + aa1 + '\t\t' + chaincode1 + '\t\t'+ atom2 + '\t\t' + aa2 + '\t\t' + chaincode2 +'\t\t' + str(distance)+ '\n')
            #DISSULFIDE_BOND----DISSULFIDE_BOND----DISSULFIDE_BOND----DISSULFIDE_BOND----DISSULFIDE_BOND----DISSULFIDE_BOND----DISSULFIDE_BOND----
            if atom1 in ligdb[:] and atom2 in ligdb[:] and 0.0 < distance < 2.2:
                db+=1
                f.write('Dissulfide_bond' + '\t\t' + atom1 + '\t\t' + aa1 + '\t\t' + chaincode1 + '\t\t'+ atom2 + '\t\t' + aa2 + '\t\t' + chaincode2 +'\t\t' + str(distance)+ '\n')
            #VAN_DER_WAALS---------VAN_DER_WAALS---------VAN_DER_WAALS---------VAN_DER_WAALS---------VAN_DER_WAALS---------VAN_DER_WAALS---------
            if aa1 in aavdw[:] and atom1 in ligvdw[:] and aa2 in aavdw[:] and atom2 in ligvdw[:] and 0.0 < distance < 5.5:
                vdw+=1
                print('{}'.format(distance))
                f.write('Van_der_Waals' + '\t\t' + atom1 + '\t\t' + aa1 + '\t\t' + chaincode1 + '\t\t'+ atom2 + '\t\t' + aa2 + '\t\t' + chaincode2 +'\t\t' + str(distance)+ '\n')
            #PI_STACKING---------PI_STACKING---------PI_STACKING---------PI_STACKING---------PI_STACKING---------PI_STACKING---------
            if aa1 in aapi[:] and aa2 in aapi[:]:
                if(chaincode1 not in Invalids and chaincode2 not in Invalids and [chaincode1,chaincode2] not in Exclusions and [chaincode2,chaincode1] not in Exclusions ):
                    coordinates1=AromaticArray[chaincode1]
                    coordinates2=AromaticArray[chaincode2]
                    aromaticdistance=np.linalg.norm(coordinates1-coordinates2)
                    if(aromaticdistance<7.2):
                        f.write('Pi_stacking  ' + '\t\t' + atom1 + '\t\t' + aa1 + '\t\t' + chaincode1 + '\t\t'+ atom2 + '\t\t' + aa2 + '\t\t' + chaincode2 +'\t\t' + str(aromaticdistance)+ '\n')
                        NormalVector1=AromaticNormals[chaincode1]/np.linalg.norm(AromaticNormals[chaincode1])
                        NormalVector2=AromaticNormals[chaincode2]/np.linalg.norm(AromaticNormals[chaincode2])
                        Angle=np.arccos(np.clip(np.dot(NormalVector1, NormalVector2), -1.0, 1.0))
                        if(Angle>50):
                            tshaped+=1
                        elif(30<Angle<50):
                            inter+=1
                        elif(Angle<30):
                            paralel+=1
                        #print(Angle,type(AromaticNormals[chaincode1]),NormalVector1,NormalVector2,np.dot(NormalVector1, NormalVector2))
                        lpi+=1
                        Exclusions.append([chaincode1,chaincode2])
            #CATION_ARYL---------#CATION_ARYL---------#CATION_ARYL---------#CATION_ARYL---------#CATION_ARYL---------
            if aa1 in aactn[:] and atom2 in ligctn[:]:
                if(chaincode1 not in Invalids and chaincode2 not in Invalids and[chaincode1,chaincode2] not in Exclusions and [chaincode2,chaincode1] not in Exclusions ):
                    coordinates1=AromaticArray[chaincode1]
                    coordinates2=np.array([New['X'].iloc[j],New['Y'].iloc[j],New['Z'].iloc[j]])
                    aromaticdistance=np.linalg.norm(coordinates1-coordinates2)
                    if(3.4 < aromaticdistance < 4.0):
                        ctn+=1
                        f.write('Cation_Aryl' + '\t\t' + 'centroid' + '\t\t' + aa1 + '\t\t' + chaincode1 + '\t\t'+ atom2 + '\t\t' + aa2 + '\t\t' + chaincode2 +'\t\t' + str(aromaticdistance)+ '\n')
                        Exclusions.append([chaincode1,chaincode2])
            elif(atom1 in ligctn[:] and aa2 in aactn[:]):
                if(chaincode1 not in Invalids and chaincode2 not in Invalids and[chaincode1,chaincode2] not in Exclusions and [chaincode2,chaincode1] not in Exclusions ):
                    coordinates1=np.array([New['X'].iloc[i],New['Y'].iloc[i],New['Z'].iloc[i]])
                    coordinates2=AromaticArray[chaincode2]
                    aromaticdistance=np.linalg.norm(coordinates1-coordinates2)
                    if(3.4 < aromaticdistance < 4.0):
                        ctn+=1
                        f.write('Cation_Aryl' + '\t\t' + atom1 + '\t\t' + aa1 + '\t\t' + chaincode1 + '\t\t'+ 'centroid' + '\t\t' + aa2 + '\t\t' + chaincode2 +'\t\t' + str(aromaticdistance)+ '\n')
                        Exclusions.append([chaincode1,chaincode2])
            #SULFUR_ARYL-----------SULFUR_ARYL-----------SULFUR_ARYL-----------SULFUR_ARYL-----------SULFUR_ARYL-----------
            if (aa1 in aaspi[:] and atom2 in ligspi1[:]):
                if(chaincode1 not in Invalids and chaincode2 not in Invalids and[chaincode1,chaincode2] not in Exclusions and [chaincode2,chaincode1] not in Exclusions ):
                    coordinates1=AromaticArray[chaincode1]
                    coordinates2=np.array([New['X'].iloc[j],New['Y'].iloc[j],New['Z'].iloc[j]])
                    aromaticdistance=np.linalg.norm(coordinates1-coordinates2)
                    if(0 < aromaticdistance < 5.3):
                        spi+=1
                        f.write('Sulfur_Aryl  ' + '\t\t' + 'centroid' + '\t\t' + aa1 + '\t\t' + chaincode1 + '\t\t'+ atom2 + '\t\t' + aa2 + '\t\t' + chaincode2 +'\t\t' + str(aromaticdistance)+ '\n')
                        Exclusions.append([chaincode1,chaincode2])
            elif(atom1 in ligspi1[:] and aa2 in aaspi[:]):
                if(chaincode1 not in Invalids and chaincode2 not in Invalids and[chaincode1,chaincode2] not in Exclusions and [chaincode2,chaincode1] not in Exclusions ):
                    coordinates1=np.array([New['X'].iloc[i],New['Y'].iloc[i],New['Z'].iloc[i]])
                    coordinates2=AromaticArray[chaincode2]
                    aromaticdistance=np.linalg.norm(coordinates1-coordinates2)
                    if(0 < aromaticdistance < 5.3):
                        spi+=1
                        f.write('Sulfur_Aryl' + '\t\t' + atom1 + '\t\t' + aa1 + '\t\t' + chaincode1 + '\t\t'+ 'centroid' + '\t\t' + aa2 + '\t\t' + chaincode2 + '\t\t' + str(aromaticdistance)+ '\n')
                        Exclusions.append([chaincode1,chaincode2])
            #AN-------------AN-------------AN-------------AN-------------AN-------------AN-------------
            if aa1 in aaan[:] and atom2 in ligan1[:]:
                if(chaincode1 not in Invalids and chaincode2 not in Invalids and[chaincode1,chaincode2] not in Exclusions and [chaincode2,chaincode1] not in Exclusions ):
                    coordinates1=AromaticArray[chaincode1]
                    coordinates2=np.array([New['X'].iloc[j],New['Y'].iloc[j],New['Z'].iloc[j]])
                    aromaticdistance=np.linalg.norm(coordinates1-coordinates2)
                    if(2.0 < aromaticdistance < 3.0):
                        an+=1
                        f.write('Anion_Aryl' + '\t\t' + 'centroid' + '\t\t' + aa1 + '\t\t' + chaincode1 + '\t\t'+ atom2 + '\t\t' + aa2 + '\t\t' + chaincode2 + '\t\t' + str(aromaticdistance)+ '\n')
                        Exclusions.append([chaincode1,chaincode2])
            elif(atom1 in ligan1[:] and aa2 in aaan[:]):
                if(chaincode1 not in Invalids and chaincode2 not in Invalids and[chaincode1,chaincode2] not in Exclusions and [chaincode2,chaincode1] not in Exclusions ):
                    coordinates1=np.array([New['X'].iloc[i],New['Y'].iloc[i],New['Z'].iloc[i]])
                    coordinates2=AromaticArray[chaincode2]
                    aromaticdistance=np.linalg.norm(coordinates1-coordinates2)
                    if(2.0 < aromaticdistance < 3.0):
                        an+=1
                        f.write('Anion_Aryl' + '\t\t' + atom1 + '\t\t' + aa1 + '\t\t' + chaincode1 + '\t\t'+ 'centroid' + '\t\t' + aa2 + '\t\t' + chaincode2 + '\t\t' + str(aromaticdistance)+ '\n')
                        Exclusions.append([chaincode1,chaincode2])

    f.close()
    f= open('output/'+os.path.split(pdbfile)[1] +".txt","w+")
    #f.write(pdbfile+";hb="+str(hb)+";sb="+str(sb)+";db="+str(db)+";lpi="+str(lpi)+";vdw="+str(vdw)+";ctn="+str(ctn)+";an="+str(an)+";spi="+str(spi))
    f.write(pdbfile+";hb={};sb={};db={};lpi={},tshaped={},inter={},paralel={};vdw={};ctn={};an={};spi={}".format(hb,sb,db,lpi,tshaped,inter,paralel,vdw,ctn,an,spi))
    f.close()
    print("---%s seconds ---" % (time.time() - start_time))
    print(pdbfile+";hb={};sb={};db={};lpi={},tshaped={},inter={},paralel={};vdw={};ctn={};an={};spi={}".format(hb,sb,db,lpi,tshaped,inter,paralel,vdw,ctn,an,spi))
