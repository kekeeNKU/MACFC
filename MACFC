from numpy import *


def auc(tlofe, ne):
    lpp = 0  
    lnp = 0  
    flag = 0  
    v = 0  
    for i in range(-1, -size(tlofe) - 1, -1):
        if tlofe[i] == ne:
            if flag == 1:
                v += lnp * lpp
                flag = 0
                lpp = 0
            lnp += 1
        else:
            if flag == 0:
                flag = 1
            lpp += 1
    v += lnp * lpp
    auc = (n0 * n1 - v) / (n0 * n1)
    return auc


def feature_selection(f, c, max_fs):
    f_auc = []  
    f_no = [str(i) for i in range(shape(f)[0])]  
    f_mtf = full((shape(f)[0], shape(f)[1]), False)
    f_ne = []  
    fl = shape(f_no)[0]
    for j in range(fl):
        argfv = argsort(f[j])  
        slofe = c[argfv]  
        ne = slofe[0]
        a = auc(slofe, ne)
        if a < 0.5:
            if c[0] == c[-1]:
                a = 1 - a
                if ne == pos:
                    ne = neg
                else:
                    ne = pos
        f_auc.append(a)
        f_ne.append(ne)

        ml = 1
        mr = 1
        for i in range(1, size(slofe)):  
            if slofe[i] == slofe[0]:
                ml += 1
            else:
                break
        for i in range(-2, -size(slofe), -1):  
            if slofe[i] == slofe[-1]:
                mr += 1
            else:
                break
        mr = size(slofe) - mr

        if slofe[0] == slofe[-1]:  
            if not slofe[0] == ne:  
                ml = 0
            else:  
                mr = size(slofe)
        f_mtf[j][argfv[ml:mr]] = True  
    # print(f_auc)
    arg_auc = argsort(-array(f_auc))  
    FName = array(f_no)[arg_auc]  
    Fvalue = array(f)[arg_auc]
    Fauc = array(f_auc)[arg_auc]
    Fne = array(f_ne)[arg_auc]
    FmTF = array(f_mtf)[arg_auc]  
    # print('SORT VALUE', Fvalue)
    # print('SORT M', FmTF)
    # print('SORT NAME', FName)
    # print('SORT AUC', Fauc)

    kk = 0
    slen = 0
    Fmcount = ones((len(FmTF[0])))  
    Fmcount = Fmcount.astype(bool)
    for i in range(fl):
        if Fauc[i] < 0.5:
            kk += 1
        Fmcount &= FmTF[i]
        if True in Fmcount:
            slen += 1
    print('Totally ', kk, ' features with auc under 0.5')

    for i in range(fl):
        if Fauc[i] < 0.5:
            continue  
        for j in range(i + 1, fl):
            if Fauc[j] < 0.5:
                continue
            nflg = 0
            if not ((FmTF[i] & FmTF[j]) == FmTF[i]).all():
                if not ((FmTF[i] & FmTF[j]) == FmTF[j]).all():
                    nflg = 1
            if nflg == 0:  
                if FmTF[i].sum() <= FmTF[j].sum():
                    Fauc[j] = -2
                else:
                    Fauc[i] = -2
                    break
    ii = []
    for i in range(fl):
        if Fauc[i] > 0.5:
            ii.append(i)
    print('Totally ' + str(fl - len(ii)) + ' features are covered and removed.')

    FName = FName[ii]
    Fvalue = Fvalue[ii]
    Fauc = Fauc[ii]
    Fne = Fne[ii]
    FmTF = FmTF[ii]

    ## start selection
    # initial
    rnk = 2
    mv_auc = Fauc[0]
    fs = [FName[0]] 
    cpms = FmTF[0]  
    fl = shape(FName)[0]

    while rnk < max_fs and mv_auc != 1:
        ft = 0
        for j in range(fl):
            if FName[j] not in fs:
                tmpFmTF = cpms & FmTF[j]
                if not ((FmTF[j] & cpms) == cpms).all():
                    mauc = 0
                    for g in fs+[FName[j]]:
                        fval = Fvalue[argwhere(FName == g)[0][0]][tmpFmTF]
                        stwlofe = array(c)[tmpFmTF]
                        argfv = argsort(fval)
                        slofe = stwlofe[argfv]
                        tauc = auc(slofe, Fne[argwhere(FName == g)[0][0]])
                        mauc += tauc
                    tmpauc = mauc / rnk
                    if tmpauc > mv_auc:
                        mv_auc = tmpauc
                        ft = j

        if mv_auc == -2 or ft == 0:
            break
        fs.append(FName[ft])
        cpms = cpms & FmTF[ft] 
        rnk += 1
        print('\nRank-' + str(rnk-1) + ' mvAUC: ' + str(mv_auc) + '  Feature set:', fs)
    return fs, mv_auc


if __name__ == '__main__':
    f = array([[1, 0.89, 0.84, 0.58, 0.5, 0.42, 0.35, 0.45, 0.36, 0.3, 0.15, 0.2, 0.25, 0.32, 0.4, 0.48, 0.53, 0.65, 0.72, 0.78],
               [0.76, 0.72, 0.7, 0.57, 0.52, 0.68, 0.63, 0.2, 0.26, 0.1, 0.71, 0.73, 0.92, 0.86, 0.55, 0.4, 0.82, 0.49, 0.33, 0.36],
               [0.67, 0.71, 0.46, 0.20, 0.16, 0.25, 0.12, 0.5, 0.38, 0.59, 0.41, 0.55, 0.33, 0.63, 0.76, 0.88, 0.29, 0.93, 0.8, 0.85]])
    c = concatenate((zeros(10), ones(10)), axis=0)
    pos, neg = set(c)
    n0, n1 = list(c).count(pos), list(c).count(neg)

    fs, mvauc = feature_selection(f, c, shape(f)[0]+1)
    


