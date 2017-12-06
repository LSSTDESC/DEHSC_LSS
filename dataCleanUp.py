import numpy as np
def dataCleanUp(dataDF):
    # input: HSC data as a pandas dataframe.
    # Two-step cleanup:
    # 1) use <>_isnull columns:  drop objects with any null column entries: the sql query apparently adds a column for many columns,
    # encoding whether the entry in the column is null or not. Then drop the null columns.
    # 2) remove the objs. with any null entries.
    print 'Given size of the dataframe: ', np.shape(dataDF)
    isNulls= [k for k in dataDF.keys() if k.__contains__('isnull')] # all the columns that signify nulls

    indToDrop= []
    for col in isNulls:
        nullInd= np.where(dataDF[col]==True)[0]
        indToDrop+=nullInd

    # find unique indexes
    indToDrop= list(set(indToDrop))

    if len(indToDrop)>0: # remove the objects with any null entries
        print 'Dropped %s rows based on isnull columns.'%len(indToDrop)
        dataDF= dataDF.drop(indToDrop)
    # drop the columns now. irrelevant.
    dataDF= dataDF.drop(isNulls, axis= 1)
    print 'Dropped %s isnull columns.'%len(isNulls)
    
    # have nans in different columns. drop row if there's a nan.
    before= len(dataDF)
    dataDF= dataDF.dropna(axis= 0)
    print 'Dropped %s rows since they contained nan entries.'%(before-len(dataDF.keys()))
    print 'Final size of the dataframe: ', np.shape(dataDF)
    print ''
    return dataDF
