"""Detect peaks in data based on their amplitude and other features."""

import numpy as np

from statsmodels.formula.api import ols

__author__ = "Marcos Duarte, https://github.com/demotu/BMC"
__version__ = "1.0.4"
__license__ = "MIT"

'''TESTING
x=np.asarray(Augerfile['S7D72'])
'''

plt.plot()

big=['O']
major=['C','Si']

chargerange=[-15,100] # range to search for largest peaks
plotcol='S7D72'
thresh=
Augerslice.plot(x='Energy',y=Augerfile.columns[1])
Augerfile.plot(x='Energy',y=Augerfile.columns[1])

chargeshift=140

# quick plotting of pair of 1D arrays (energy vs counts or deriv)
thisen=energy[lowind[0]:highind[1]+1] # x values
plt.plot(thisen, s7d7)
plt.plot(thisen, rawdata)

def detectmultiplepeaks(Augerfile, plotcol, Elements, AESquantparams, chargerange, thresh):
    ''' Look for major peaks in smdiff data at shifted locations and compare 
    energy intervals to major peaks  '''
    # Peak search around 
    foundpeaks={} # elem, ev, val
    for elem in Elements:
        match=AESquantparams[AESquantparams['element']==elem]
        if len(match)==1:
            negpeak=match.iloc[0]['negpeak']
            searchwidth=int(match.iloc[0]['searchwidth'])
            pospeak=match.iloc[0]['pospeak']
        Augerslice=Augerfile[( Augerfile['Energy'] > negpeak-searchwidth+chargerange[0]) & ( Augerfile['Energy'] < negpeak+chargerange[1])]
        # just look at min
        evals=detect_peaks(Augerslice, mph=None, mpd=1, threshold=0, edge='rising',
                 kpsh=False, valley=True, show=True, ax=None)
        # Peak containing subset of df
        regpeaks=Augerslice[Augerslice['Energy'].isin(evals)]
        


def findvalidpeaks(regpeaks, thresh):
    ''' Looking for peaks significantly large than background '''
    # outlier test on deriv peaks
    xdata=regpeaks['Energy'].tolist()
    ydata=regpeaks[regpeaks.columns[1]].tolist()
    regression= ols("data ~ x", data=dict(data=ydata, x=xdata)).fit()
    outliers=regression.outlier_test() 
    colnames=['Resid','Pval','Bonf']
    outliers.columns=colnames
        
df=Augerslice    
    
def findenergies(Elements, AESquantparams):
    '''  '''

def detect_peaks(df, mph=None, mpd=1, threshold=0, edge='rising',
                 kpsh=False, valley=True, show=True, ax=None):

    """Detect peaks in data based on their amplitude and other features.

    Parameters
    ----------
    x : 1D array_like
        data.
    mph : {None, number}, optional (default = None)
        detect peaks that are greater than minimum peak height (just y value)
    mpd : positive integer, optional (default = 1)
        detect peaks that are at least separated by minimum peak distance (in
        number of data).
    threshold : positive number, optional (default = 0)
        detect peaks (valleys) that are greater (smaller) than `threshold`
        in relation to their immediate neighbors.
    edge : {None, 'rising', 'falling', 'both'}, optional (default = 'rising')
        for a flat peak, keep only the rising edge ('rising'), only the
        falling edge ('falling'), both edges ('both'), or don't detect a
        flat peak (None).
    kpsh : bool, optional (default = False)
        keep peaks with same height even if they are closer than `mpd`.
    valley : bool, optional (default = False)
        if True (1), detect valleys (local minima) instead of peaks.
    show : bool, optional (default = False)
        if True (1), plot data in matplotlib figure.
    ax : a matplotlib.axes.Axes instance, optional (default = None).
 
    Returns
    -------
    ind : 1D array_like
        indeces of the peaks in `x`.

    Notes
    -----
    The detection of valleys instead of peaks is performed internally by simply
    negating the data: `ind_valleys = detect_peaks(-x)`
    
    The function can handle NaN's 

    See this IPython Notebook [1]_.

    References
    ----------
    .. [1] http://nbviewer.ipython.org/github/demotu/BMC/blob/master/notebooks/DetectPeaks.ipynb

    Examples
    --------
    >>> from detect_peaks import detect_peaks
    >>> x = np.random.randn(100)
    >>> x[60:81] = np.nan
    >>> # detect all peaks and plot data
    >>> ind = detect_peaks(x, show=True)
    >>> print(ind)

    >>> x = np.sin(2*np.pi*5*np.linspace(0, 1, 200)) + np.random.randn(200)/5
    >>> # set minimum peak height = 0 and minimum peak distance = 20
    >>> detect_peaks(x, mph=0, mpd=20, show=True)

    >>> x = [0, 1, 0, 2, 0, 3, 0, 2, 0, 1, 0]
    >>> # set minimum peak distance = 2
    >>> detect_peaks(x, mpd=2, show=True)

    >>> x = np.sin(2*np.pi*5*np.linspace(0, 1, 200)) + np.random.randn(200)/5
    >>> # detection of valleys instead of peaks
    >>> detect_peaks(x, mph=0, mpd=20, valley=True, show=True)

    >>> x = [0, 1, 1, 0, 1, 1, 0]
    >>> # detect both edges
    >>> detect_peaks(x, edge='both', show=True)

    >>> x = [-2, 1, -2, 2, 1, 1, 3, 0]
    >>> # set threshold = 2
    >>> detect_peaks(x, threshold = 2, show=True)
    """
    xvals=df.Energy
    xvals=xvals.reset_index(drop=True)
    yvals=np.asarray(df[df.columns[1]])

    yvals = np.atleast_1d(yvals).astype('float64')
    if yvals.size < 3: # length of data
        return np.array([], dtype=int) # return empty list if too short
    if valley: # flips data
        yvals = -yvals
    # find indices of all peaks
    dx = yvals[1:] - yvals[:-1]
    # handle NaN's
    indnan = np.where(np.isnan(yvals))[0]
    if indnan.size: # sets nan values to 
        yvals[indnan] = np.inf
        dx[np.where(np.isnan(dx))[0]] = np.inf
    ine, ire, ife = np.array([[], [], []], dtype=int)
    if not edge:  # handles flat peaks
        ine = np.where((np.hstack((dx, 0)) < 0) & (np.hstack((0, dx)) > 0))[0]
    else:
        if edge.lower() in ['rising', 'both']:
            ire = np.where((np.hstack((dx, 0)) <= 0) & (np.hstack((0, dx)) > 0))[0]
        if edge.lower() in ['falling', 'both']:
            ife = np.where((np.hstack((dx, 0)) < 0) & (np.hstack((0, dx)) >= 0))[0]
    # list of indices with detected peak
    ind = np.unique(np.hstack((ine, ire, ife)))
    # handle NaN's
    if ind.size and indnan.size:
        # NaN's and values close to NaN's cannot be peaks
        ind = ind[np.in1d(ind, np.unique(np.hstack((indnan, indnan-1, indnan+1))), invert=True)]
    # first and last values of yvals cannot be peaks (removes these indices)
    if ind.size and ind[0] == 0:
        ind = ind[1:]
    if ind.size and ind[-1] == yvals.size-1:
        ind = ind[:-1]
    # remove peaks < minimum peak height
    if ind.size and mph is not None:
        ind = ind[yvals[ind] >= mph]
    # remove peaks - neighbors < threshold
    if ind.size and threshold > 0:
        dx = np.min(np.vstack([yvals[ind]-yvals[ind-1], yvals[ind]-yvals[ind+1]]), axis=0)
        ind = np.delete(ind, np.where(dx < threshold)[0])
    # detect small peaks closer than minimum peak distance
    if ind.size and mpd > 1:
        ind = ind[np.argsort(yvals[ind])][::-1]  # sort ind by peak height
        idel = np.zeros(ind.size, dtype=bool)
        for i in range(ind.size):
            if not idel[i]:
                # keep peaks with the same height if kpsh is True
                idel = idel | (ind >= ind[i] - mpd) & (ind <= ind[i] + mpd) \
                    & (yvals[ind[i]] > yvals[ind] if kpsh else True)
                idel[i] = 0  # Keep current peak
        # remove the small peaks and sort back the indices by their occurrence
        ind = np.sort(ind[~idel])
    evals=xvals[xvals.index.isin(ind.tolist())]
    evals=[int(i) for i in evals]
    if show:
        if indnan.size:
            yvals[indnan] = np.nan
        if valley:
            yvals = -yvals
        _plot(yvals, mph, mpd, threshold, edge, valley, ax, ind)

    return evals


def _plot(yvals, mph, mpd, threshold, edge, valley, ax, ind):
    """Plot results of the detect_peaks function, see its help."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print('matplotlib is not available.')
    else:
        if ax is None:
            _, ax = plt.subplots(1, 1, figsize=(8, 4))

        ax.plot(yvals, 'b', lw=1)
        if ind.size:
            label = 'valley' if valley else 'peak'
            label = label + 's' if ind.size > 1 else label
            ax.plot(ind, yvals[ind], '+', mfc=None, mec='r', mew=2, ms=8,
                    label='%d %s' % (ind.size, label))
            ax.legend(loc='best', framealpha=.5, numpoints=1)
        ax.set_xlim(-.02*yvals.size, yvals.size*1.02-1)
        ymin, ymax = yvals[np.isfinite(yvals)].min(), yvals[np.isfinite(yvals)].max()
        yrange = ymax - ymin if ymax > ymin else 1
        ax.set_ylim(ymin - 0.1*yrange, ymax + 0.1*yrange)
        ax.set_xlabel('Data #', fontsize=14)
        ax.set_ylabel('Amplitude', fontsize=14)
        mode = 'Valley detection' if valley else 'Peak detection'
        ax.set_title("%s (mph=%s, mpd=%d, threshold=%s, edge='%s')"
                     % (mode, str(mph), mpd, str(threshold), edge))
        # plt.grid()
        plt.show()