from numpy import *
import data_process as dp


def auc(tlofe, ne):
    lpp = 0
    lnp = 0
    flag = 0
    aac = 0
    for i in range(-1, -size(tlofe) - 1, -1):
        if tlofe[i] == ne:
            if flag == 1:
                aac += lnp * lpp
                flag = 0
                lpp = 0
            lnp += 1
        else:
            if flag == 0:
                flag = 1
            lpp += 1
    aac += lnp * lpp
    auc = (n0 * n1 - aac) / (n0 * n1)
    return auc


def feature_selection(f, c, max_fs):
    f_auc = []
    f_no = [str(i) for i in range(shape(f)[0])]
    f_mtf = full((shape(f)[0], shape(f)[1]), False)
    f_ne = []
    fl = shape(f_no)[0]
    # print('fl ', fl)
    for j in range(fl):
        argfv = argsort(f[j])
        slofe = c[argfv]
        ne = slofe[0]
        a = auc(slofe, ne)
        if a < 0.5:
            if slofe[0] == slofe[-1]:
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
    # print('Totally ', kk, ' features with auc under 0.5')

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
    # print('Totally ' + str(fl - len(ii)) + ' features are covered and removed.')

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
        temp = 0
        for j in range(fl):
            if FName[j] not in fs:
                tmpFmTF = cpms & FmTF[j]
                if not ((FmTF[j] & cpms) == cpms).all():
                    mauc = 0
                    for g in fs + [FName[j]]:
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
                        temp = Fauc[j]
                    elif tmpauc == mv_auc and Fauc[j] > temp:
                        ft = j
                        temp = Fauc[j]


        if mv_auc == -2 or ft == 0:
            break
        fs.append(FName[ft])
        cpms = cpms & FmTF[ft]
        rnk += 1
        print('\nRank-' + str(rnk - 1) + ' mvAUC: ' + str(mv_auc) + '  Feature set:', fs)
    return fs, mv_auc



if __name__ == '__main__':
    # load UCI data
    f, c = dp.load_ult()

    # load TCGA data
    # fname, f, c = dp.load_TCGA('./TCGA/genomicMatrix-PRAD.xls')

    pos, neg = set(c)
    n0, n1 = list(c).count(pos), list(c).count(neg)

    fs, mvauc = feature_selection(f, c, shape(f)[0] + 1)
    # print(fname[array(fs).astype(int)])   # print gene names for TCGA data



