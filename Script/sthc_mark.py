from scipy.spatial.distance import pdist, squareform
from scipy.spatial import distance
from scipy import spatial
from scipy.cluster.hierarchy import dendrogram
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import fcluster
from sklearn.metrics import silhouette_samples, silhouette_score
from shapely.geometry.point import Point
from datetime import timedelta
from dateutil import parser
import datetime
import pandas as pd
import numpy as np
from pyproj import Transformer
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy.cluster import hierarchy
from scipy import spatial
from scipy.ndimage import gaussian_filter
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import pdist


class MissingDataset(Exception):
    """Exception raised for errors in the input.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message


class STHC(object):
    def __init__(self):
        self.functions = []
        self.data = None
        self.df = None
        self.X = None
        self.distance_params = {}
        self.timedistance_params = {}
        self.attribute_params = {}
        self.tree = None
        self.attribute_norm_values = {}
        self.sqdm = None

    def calculate_clusters(self,method='average'):
        self.row_clusters = linkage(self.dm, method=method)
        c,co_dists = cophenet(self.row_clusters,self.dm)
        print (method+' Cophenetic correlation coeffficient = '+str(c))


    def print_dendrogram(self, figsize=(9, 6), dpi=300):
        colors_grey = [n for n, hex in colors.cnames.iteritems() if 'grey' in n]
        # print colors_grey
        hierarchy.set_link_color_palette(None)  # hierarchy.set_link_color_palette(colors_grey)
        color_list = plt.get_cmap('Greys')
        fig = plt.figure(figsize=figsize, dpi=dpi)

        # frame1.axes.yaxis.set_ticklabels([])
        plt.title('Hierarchical Clustering Dendrogram')
        plt.xlabel('Index')
        plt.ylabel('Distance')
        row_dendr = dendrogram(self.row_clusters)
        frame1 = plt.gca()
        # frame1.set_position([0, 0, 1, 1])
        frame1.axes.xaxis.set_ticklabels([])
        plt.show()

    def calculate_distance_matrix(self):

        if isinstance(self.df, type(None)):
            raise MissingDataset("calcualte_distance_matrix", "Need to load the dataset first...")
        else:
            # self.dm_old = pdist(self.X, self.calculate_distance)
            try:
                del self.dm
            except:
                pass
            self.dm = self.pdist_custom(self.X, self.distfunc, self.timefunc, self.attribute_params['function'],
                                        self.weights)

    def calculate_distance(self, x, y):
        weights = []
        finalsum = []
        # Calculate distance between two points

        dist = self.distfunc(x[0:2], y[0:2])
        dist_norm = dist / x[2]
        # print(dist_norm)
        finalsum.append(dist_norm)
        weights.append(self.distance_params['weights'])
        # calcualte the difference in time

        dist_time = self.timefunc(x[3], y[3])
        dist_time_norm = dist_time / x[4]
        # print(dist_time_norm)

        finalsum.append(dist_time_norm)
        weights.append(self.timedistance_params['weights'])
        # Calculate values for attributes
        # dist_att = []
        dist_att_norm = []
        for i in range(0, len(self.attribute_params['fieldnames']) + 2, 2):
            x_i = i + 5
            x_n_i = i + 5 + 1
            dist_att_norm.append((x[x_i], y[x_i], x[x_n_i]))
        for i in range(0, len(self.attribute_params['fieldnames'])):
            dist_a_i = self.attribute_params['function'][i](dist_att_norm[i][0], dist_att_norm[i][1])
            # print(dist_a_i)
            finalsum.append(dist_a_i / dist_att_norm[i][2])
            weights.append(self.attribute_params['weights'][i])

        # print(finalsum)

        res = np.sum([v * w for v, w in zip(finalsum, weights)]) / np.sum(weights)

        if res >= 0.0:
            return res
        else:
            return 0.0
        # return 1.0

    def pdist_custom(self, X, distfunc, timefunc, attributefuncs, weights):
        """distance matrix that doesn't rely on numpy's pdist. speed check"""
        n = X.shape[0]

        # empty array
        dm = np.zeros(int(n * (n - 1) / 2.0))
        for i in range(n - 1):
            x = X[i]
            y = X[i + 1:]
            ress = []
            # for y in Y:
            finalsum = 0.0

            if distfunc:
                dist = distfunc(x[0:2], y[:, 0:2])[0]
                # print(dist)
                dist_norm = dist / x[2]
                # print(dist_norm)
                finalsum = dist_norm * weights[0]

            if timefunc:
                dist_time = timefunc(x[3], y[:, 3])
                # print(dist_time)
                dist_time_norm = dist_time / x[4]
                # print(dist_time_norm)
                finalsum += dist_time_norm * weights[1]

            base = 5
            for j in range(0, len(attributefuncs)):
                x_i = j + base
                x_n_i = j + base + 1
                f = attributefuncs[j]
                if f != None:
                    t = f(x[x_i], y[:, x_i]) / x[x_n_i] * weights[j + 2]
                    finalsum = finalsum + t
                base += 1
            res = finalsum / np.sum(weights)
            # ress.append(res)
            idx_s = int((2 * n - i - 1) * i / 2.0)
            idx_e = int((2 * n - i - 2) * (i + 1) / 2.0)
            dm[idx_s:idx_e] = res

        return dm

    def load_dataset_csv(self, path, params, toepsg=None, dayfirst=True, sample=None):
        """Load a csv file to be used in the space-time hierarchical clustering
        path: string path to the csv file to be used.
        params: dictionary specifying the functions to calculate the distance between different components
        toepsg: epsg code to transform latitude and longitude coordinates to planar coordinates.
        dayfirst: for date conversion. if days are first in the date string.
        sample: if None, use the full dataset, if <1 use a percentage of the dataset, if > 1 use as the sample size
        Planar coordinates are used because the points are placed in a kdtree for faster calculations.
        """
        print("Reading the file...")
        self.df = pd.read_csv(path)
        if sample:
            if sample < 1:
                self.df = self.df.sample(frac=sample, replace=False)
            elif sample > 1:
                self.df = self.df.sample(n=sample, replace=False)
        self.distance_params = params['distance'].copy()
        self.timedistance_params = params['timedistance'].copy()
        self.attribute_params = params['attributes'].copy()
        self.distance_params['field_index'] = [0,
                                               1]  # need to now the structure of the array when calculating between points
        self.timedistance_params['field_index'] = [2]
        # self.attribute_params['field_index'] = [i+3 for i,x in enumerate(self.attribute_params['fieldnames'])]
        self.distfunc = self.distance_params['function']
        self.timefunc = self.timedistance_params['function']
        fieldlist = []  # for extracting the fields into the array used in the distance calculations

        xfn = self.distance_params['fieldnames']['x']
        yfn = self.distance_params['fieldnames']['y']
        # dtfn = self.timedistance_params['fieldnames']['datetime']
        dtfn = self.timedistance_params['fieldnames']['datetime']
        attributefields = self.attribute_params['fieldnames']

        xs = []  # if converting over
        ys = []
        dts = []
        print(xfn)
        print(yfn)
        print('Processing table for datetime and coordinates...')

        # for idx,row in self.df.iterrows():
        #    if toepsg:
        #        try:
        #            xt, yt = self.transform_pnt(toepsg,row[xfn],row[yfn])
        #        except:
        #            xt = np.nan
        #            yt = np.nan
        #        xs.append(xt)
        #        ys.append(yt)
        #    try:
        #        dtt = self.convertDate(row[dtfn])
        #    except:
        #        dtt = None

        #    dts.append(dtt)
        #    if idx%1000==0:
        #        print(idx)

        if toepsg:
            transformer = Transformer.from_crs(4326, toepsg)
            res = self.df[[xfn, yfn]].apply(lambda x: transformer.transform(x[0], x[1]), axis=1)
            self.df[xfn] = res.apply(lambda x: x[0])
            self.df[yfn] = res.apply(lambda x: x[1])

        if self.distance_params['normalize_n']:
            self.df['dist_norm'] = np.ones(len(self.df)) * self.distance_params['normalize_n']
        elif self.distance_params['normalize_k']:
            self.tree = spatial.KDTree(self.df[[xfn, yfn]].values)
            nn = self.df[[xfn, yfn]].apply(lambda x: self.tree.query(x, k=self.distance_params['normalize_k']), axis=1)
            self.df['dist_norm'] = nn.apply(lambda x: np.max(x[0]))
            self.df['dist_norm'] = self.df['dist_norm'] + 1.0
        elif self.distance_params['normalize_mean']:
            raise NotImplementedError(
                "normalizing distannce by mean of distance not implemented, use normalize_n or nearest neighbor")

        # Old time distance implementation -> need to convert if YearQuarter is used
        # self.df[dtfn] = self.df[dtfn].apply(lambda x: parser.parse(x,dayfirst=dayfirst))

        # if self.timedistance_params['normalize_n']:
        #    self.df['timenorm'] =  np.ones(len(self.df))*self.timedistance_params['normalize_n']
        # elif self.timedistance_params['normalize_k']:
        #    raise NotImplementedError("normalizing by k nearest neighbors for time not implemented, use normalize_n")
        # elif  self.distance_params['normalize_mean'] == True:
        #    raise NotImplementedError("normalizing by mean of time not implemented, use normalize_n")

        # YearQuarter timedistance
        if self.timedistance_params['normalize_n']:
            self.df['timenorm'] = np.ones(len(self.df)) * self.distance_params['normalize_n']
        elif self.timedistance_params['normalize_k']:
            self.tree = spatial.KDTree(self.df[[xfn, yfn]].values)
            nn = self.df[[xfn, yfn]].apply(lambda x: self.tree.query(x, k=self.distance_params['normalize_k']), axis=1)
            self.df['timenorm'] = nn.apply(lambda x: np.max(x[0]))
            self.df['timenorm'] = self.df['dist_norm'] + 1.0
        elif self.timedistance_params['normalize_mean']:
            raise NotImplementedError(
                "normalizing distannce by mean of distance not implemented, use normalize_n or nearest neighbor")

        fieldlist.append(xfn)
        fieldlist.append(yfn)
        fieldlist.append('dist_norm')
        fieldlist.append(dtfn)
        fieldlist.append('timenorm')

        print('Preparing attribute normalization....')
        normalize_fields = []
        for i, afn in enumerate(self.attribute_params['fieldnames']):
            if self.attribute_params['normalize_n'][i]:
                self.df['attnorm{0}'.format(i)] = np.ones(len(self.df)) * self.attribute_params['normalize_n'][i]
                normalize_fields.append('attnorm{0}'.format(i))
            elif self.attribute_params['normalize_mean'][i] == True:
                self.df['attnorm{0}'.format(i)] = np.ones(len(self.df)) * self.df[afn].mean()
                normalize_fields.append('attnorm{0}'.format(i))
            elif self.attribute_params['normalize_k']:
                raise NotImplementedError(
                    "normalizing by k nearest neighbors for attributes not implemented, use normalize_n")

        for af, afn in zip(attributefields, normalize_fields):
            fieldlist.append(af)
            fieldlist.append(afn)

        self.X = self.df[fieldlist].values

        self.weights = []
        self.weights.append(params['distance']['weights'])
        self.weights.append(params['timedistance']['weights'])
        self.weights += params['attributes']['weights']

    def elbow_plot(self, n_clusters=None):
        if n_clusters == None:
            n_clusters = list(range(2, int(len(self.X) / 2), 10))
        outScores = []
        for n_c in n_clusters:
            try:
                # print "number of clusters %s"%n_c
                cluster_labels = fcluster(self.row_clusters, n_c, criterion='maxclust')
                # print cluster_labels
                ssq = self.silhouette_avg_custom(self.X.shape[0], self.dm, cluster_labels, self.square_to_condensed)
                # silhouette_avg = silhouette_score(self.sqdm + self.sqdm.T,cluster_labels,metric='euclidean')
                outScores.append(ssq)
            except:
                pass
        a = np.array(list(zip(n_clusters, outScores)))
        smoothed = gaussian_filter(a[:, 1], 3.)
        bottom = np.gradient(smoothed)[::-1].argmin()
        # print(bottom)
        # print(len(smoothed) - np.gradient(smoothed)[::-1].argmin())
        fig = plt.figure(figsize=(5, 5), dpi=300)
        ax = fig.add_subplot(111)
        ax.plot(a[:, 0], a[:, 1], '.')
        plt.plot(a[:, 0], smoothed, '-')
        plt.plot(a[bottom][0], smoothed[bottom], 'o')
        ax.set_xticks(n_clusters)
        print("Try {0} clusters...".format(a[bottom][0]))
        # print("Try {0} clusters...".format(a[len(smoothed)-np.gradient(smoothed)[::-1].argmin()][0]))
        plt.show()
        return list(zip(n_clusters, outScores))

    def elbow_ssq(self, n, dm, cluster_labels, vf):
        f = np.vectorize(vf)
        scores = 0.0
        for c in np.unique(cluster_labels):
            indexes = (cluster_labels == c).nonzero()[0]
            all_distances = []
            for i in range(len(indexes)):
                x = indexes[i]
                y = indexes[i + 1:]
                dist_within = np.array([])
                if i != 0:
                    y = np.append(y, indexes[0:i])
                if len(y) > 0:
                    dist_within = dm[f(x, y, n).astype(int)]
                all_distances += list(dist_within)
            scores += np.var(all_distances)
        return scores

    def elbow_alt(self, max_n_clusters=None):
        # from here:https://joernhees.de/blog/2015/08/26/scipy-hierarchical-clustering-and-dendrogram-tutorial/
        if max_n_clusters == None:
            max_n_clusters = int(len(self.X) / 2)
        last = self.row_clusters[-max_n_clusters:, 2]
        last_rev = last[::-1]
        idxs = np.arange(1, len(last) + 1)
        plt.plot(idxs, last_rev)

        slope = np.diff(last, 1)
        slope = np.gradient(last)
        slope_rev = slope[::-1]
        # print ("Optimal Clusters Slope: {0}".format(slope_rev.argmax()+2))
        acceleration = np.diff(last, 2)  # 2nd derivative of the distances
        acceleration_rev = acceleration[::-1]
        plt.plot(idxs, np.gradient(last)[::-1])
        k = len(last) - np.gradient(last)[::-1].argmin() + 2
        plt.plot(last[k], k, ".")
        plt.show()
        # k = acceleration_rev.argmax() + 2  # if idx 0 is the max of this we want 2 clusters
        print("Try {0} clusters...".format(len(last) - np.gradient(last)[::-1].argmin() + 2))
        return acceleration

    def square_to_condensed(self, i, j, n):
        assert i != j, "no diagonal elements in condensed matrix"
        if i < j:
            i, j = j, i
        return int(n * j - j * (j + 1) // 2 + i - 1 - j)

    def silhouette_avg_custom(self, n, dm, cluster_labels, vf):
        f = np.vectorize(vf)
        scores = []
        for c in np.unique(cluster_labels):
            indexes = (cluster_labels == c).nonzero()[0]
            for i in range(len(indexes)):
                x = indexes[i]
                y = indexes[i + 1:]
                dist_within = np.array([])
                if i != 0:
                    y = np.append(y, indexes[0:i])
                if len(y) > 0:
                    dist_within = dm[f(x, y, n).astype(int)]
                    dist_outside = np.array([])
                    for c_j in np.unique(cluster_labels[cluster_labels != c]):
                        indexes_not = (cluster_labels == c_j).nonzero()[0]
                        dist_outside = np.append(dist_outside, np.mean(dm[f(x, indexes_not, n).astype(int)]))
                    scores.append(
                        (np.min(dist_outside) - np.mean(dist_within)) / max(np.mean(dist_within), np.min(dist_outside)))
                else:
                    scores.append(0)

        return np.mean(scores)

    def silhoutte_scores_n_clusters(self, n_clusters=None, printV=False):
        outScores = []
        if isinstance(self.sqdm, type(None)):
            self.sqdm = squareform(self.dm)
        if n_clusters == None:
            n_clusters=list(range(1,50))
            #n_clusters= list(range(1,len(self.X)))
            #print(str(len(self.X))+'potential clusters')
            #n_clusters = list(range(2, int(len(self.X) / 2), 10))
        for n_c in n_clusters:
            print(str(n_c)+'processed')
            try:
                # print "number of clusters %s"%n_c
                cluster_labels = fcluster(self.row_clusters, n_c, criterion='maxclust')
                # print cluster_labels
                silhouette_avg = silhouette_score(self.sqdm, cluster_labels, sample_size = 1000)
                #silhouette_avg = self.silhouette_avg_custom(self.X.shape[0], self.dm, cluster_labels,
                 #                                               self.square_to_condensed)
                # silhouette_avg = silhouette_score(self.sqdm + self.sqdm.T,cluster_labels,metric='euclidean')
                outScores.append(silhouette_avg)
                    #if printV:
                print("For n_clusters =%s" % n_c)
                print("The average silhouette_score is : %s" % silhouette_avg)

            except:
                pass
        #meanscore = np.mean(outScores)
        # print(meanscore)
        #smoothed = #gaussian_filter(outScores, 2.)
        #bottom = len(smoothed) - np.gradient(smoothed)[::-1].argmin()
        # print(np.gradient(smoothed)[::-1].argmin())
        #fig = plt.figure(figsize=(5, 5), dpi=300)
        #ax = fig.add_subplot(111)
        #ax.plot(n_clusters, outScores, '.')
        #plt.plot(n_clusters, smoothed, '-')
        #plt.plot(n_clusters[bottom], smoothed[bottom], 'o')
        #ax.set_xticks(n_clusters)
        #print("Try {0} clusters...".format(np.array(n_clusters)[np.array(outScores) > meanscore][0]))
        #print("Try {0} clusters...".format(n_clusters[bottom]))
        #plt.show()
        return list(zip(n_clusters, outScores))

    def silouhetteScore(self, n_clusters, printV=False):
        outScores = []
        if isinstance(self.sqdm, type(None)):
            self.sqdm = squareform(self.dm)
        for n_c in n_clusters:
            try:
                # print "number of clusters %s"%n_c
                cluster_labels = fcluster(self.row_clusters, n_c, criterion='maxclust')
                # print cluster_labels
                silhouette_avg = silhouette_score(self.sqdm, cluster_labels, metric='precomputed')
                # silhouette_avg = silhouette_score(self.sqdm + self.sqdm.T,cluster_labels,metric='euclidean')
                outScores.append(silhouette_avg)
                if printV:
                    print("For n_clusters =%s" % n_c)
                    print("The average silhouette_score is : %s" % silhouette_avg)
            except:
                pass
        return outScores

    def plotClusters3D(self):
        plt.style.use('seaborn-whitegrid')
        mpl.rc('font', family='Times New Roman')
        mpl.rcParams.update({'font.size': 7})
        fig1 = plt.figure(figsize=(6, 6), dpi=300)
        ax1 = fig1.add_subplot(111, projection='3d')
        ax1.set_aspect('equal')
        ax1.view_init(20, 60)
        # gdf_eqwt.plot(ax=ax,column='c2_lbl',cmap='Dark2', categorical=True,markersize=6,marker="^")
        classes = {1: {'xs': [], 'ys': [], 'zs': [], 'c': 'k', 'm': 'o'},
                   2: {'xs': [], 'ys': [], 'zs': [], 'c': 'gray', 'm': '^'}}

        ax1.xaxis.set_major_formatter(ticker.FuncFormatter(meters_formatter))
        ax1.yaxis.set_major_formatter(ticker.FuncFormatter(meters_formatter))

        ax1.xaxis.set_major_locator(ticker.MultipleLocator(200))
        ax1.yaxis.set_major_locator(ticker.MultipleLocator(200))

        ax1.set_zticks([0, 14400, 28800, 43200, 57600, 72000, 86400])
        for idx, row in gdf_wttime.iterrows():
            classes[row['c2_lbl']]['xs'].append(row['geometry'].x)
            classes[row['c2_lbl']]['ys'].append(row['geometry'].y)
            classes[row['c2_lbl']]['zs'].append(row['Time'])

        for k, v in classes.iteritems():
            ax1.scatter(v['xs'], v['ys'], v['zs'], c=v['c'], marker=v['m'])
            # ax1.plot(v['xs'],v['ys'],v['zs'],color=colorClusters[k],alpha=.8,marker='o',linestyle='')
            # ax1.plot(v['xs'],v['ys'],v['zs'],alpha=.8,marker='o',linestyle='')
            # plt.plot()
            # plt.plot(ks,minSc)
            # plt.plot(ks,tenc)
        ax1.set_zticklabels(["0:00", "4:00", "8:00", "12:00", "16:00", "20:00", "23:59"])
        ax1.set_xlabel('X')
        ax1.set_ylabel('Y')
        ax1.set_zlabel('Time of Day')

