import numpy as np
import pandas as pd
from typing import List, Tuple
from networks import create_region_species_network
from complexity import compute_eci, compute_fitness, count
import plotly.graph_objects as go
from scipy.stats import linregress
from plotting_andrea import EcosystemComplexityData, Plotter, RegionData


"""

THIS IS SHIT CODE, DO NOT RUN IT

"""


class CleanData:
    def __init__(self):
        path = '/Users/andreamusso/Desktop/PhD/CSSS/TreeSpace/GitHub/region_data/temperatures.csv'
        df = pd.read_csv(path, header=0)
        rd = RegionData(geo_groups=['Europe'])
        list_data = []
        for region in rd.df_region_data.region.values:
            try:
                avg_temp, min_temp, max_temp = self.get_data(df=df, region=region)
                list_data.append([avg_temp, min_temp, max_temp])
            except:
                print('REGION', region)
                if region == 'Macedonia':
                    avg_temp, min_temp, max_temp = self.get_data(df=df, region='North Macedonia')
                    list_data.append([avg_temp, min_temp, max_temp])
                if region == 'Moldova':
                    avg_temp, min_temp, max_temp = self.get_data(df=df, region='Republic of Moldova')
                    list_data.append([avg_temp, min_temp, max_temp])
                if region == 'Czech Republic':
                    avg_temp, min_temp, max_temp = self.get_data(df=df, region='Czechia')
                    list_data.append([avg_temp, min_temp, max_temp])
                if region == 'San Marino':
                    avg_temp, min_temp, max_temp = self.get_data(df=df, region='Slovenia')
                    list_data.append([avg_temp, min_temp, max_temp])
                if region == 'Kosovo':
                    avg_temp, min_temp, max_temp = self.get_data(df=df, region='Albania')
                    list_data.append([avg_temp, min_temp, max_temp])



        list_data = np.array(list_data)
        rd.df_region_data['avg_temp'] = list_data[:, 0]
        rd.df_region_data['min_temp'] = list_data[:, 1]
        rd.df_region_data['max_temp'] = list_data[:, 2]
        rd.df_region_data.to_csv(path_or_buf='/Users/andreamusso/Desktop/PhD/CSSS/TreeSpace/GitHub/region_data/EuropeData.csv')

    def get_data(self, df, region):
        df_region = df.loc[df['Country'] == region]
        avg_temp = df_region['avg_temp'].values[0]
        min_temp = df_region['min_temp'].values[0]
        max_temp = df_region['max_temp'].values[0]
        return avg_temp, min_temp, max_temp


class PlotRegionComplexity:
    def __init__(self, geo_groups: List[str]):
        """
        :param geo_groups: the geographic region groups (e.g. Europe, China..) for which you want to load the data
        """

        self.region_data = RegionData(geo_groups=geo_groups)
        self.eco_complexity_data = EcosystemComplexityData(geo_groups=geo_groups)

    def plot_complexity_vs_var(self, complexity_measure: str, xvar: str, color_var: str = ''):
        df_region, df_complexity = self.get_cleaned_data()
        x = df_region[xvar].values.astype(float)
        y = df_complexity[complexity_measure].values.astype(float)
        if color_var != '':
            marker_color = df_region[color_var].values.astype(float)
        else:
            marker_color = 'Black'

        text = df_region['region'].values.astype(str)
        Plotter.plot(x=x, y=y, text=text, xlabel=xvar, ylabel=f'Complexity ({complexity_measure})', xscale='linear',
                     marker_color=marker_color, yscale='linear', title='')

    def get_cleaned_data(self):
        very_small_regions = ['Monaco', 'San Marino', 'Liechtenstein', 'Luxembourg', 'Andorra', 'Vatican City', 'Macao', 'Hong Kong', 'Shanghai']
        df_region_data = self.region_data.df_region_data[~(self.region_data.df_region_data.region.isin(very_small_regions))].copy()
        df_region_data.sort_values('region', inplace=True)
        df_complexity = self.eco_complexity_data.df_regions[~(self.eco_complexity_data.df_regions.region.isin(very_small_regions))].copy()
        df_complexity.sort_values('region', inplace=True)
        return df_region_data, df_complexity


class PlotMaps:
    def __init__(self, geo_groups):
        """
        :param geo_groups: the geographic region groups (e.g. Europe, China..) for which you want to load the data
        """

        self.region_data = RegionData(geo_groups=geo_groups)
        self.eco_complexity_data = EcosystemComplexityData(geo_groups=geo_groups)

    def plot_maps(self):
        import plotly.express as px
        import json


        df, _ = self.get_cleaned_data()
        fig = px.choropleth(locations=df.region.values, locationmode="country names", scope="europe",
                            color=df['avg_temp'].values.astype(float),
                            color_continuous_scale=px.colors.sequential.Plasma,
                            labels={'color': 'Average Temperature'}
                            )

        fig.show(renderer="browser")

    def get_cleaned_data(self):
        very_small_regions = ['Monaco', 'San Marino', 'Liechtenstein', 'Luxembourg', 'Andorra', 'Vatican City', 'Macao', 'Hong Kong', 'Shanghai']
        df_region_data = self.region_data.df_region_data[~(self.region_data.df_region_data.region.isin(very_small_regions))].copy()
        df_region_data.sort_values('region', inplace=True)
        df_complexity = self.eco_complexity_data.df_regions[~(self.eco_complexity_data.df_regions.region.isin(very_small_regions))].copy()
        df_complexity.sort_values('region', inplace=True)
        return df_region_data, df_complexity



if __name__ == '__main__':
    geo_groups = ['Europe']

    plots = PlotRegionComplexity(geo_groups=geo_groups)

    complexity_measures = ['eci_r', 'fit_r', 'cou_r']
    for cm in complexity_measures:
        plots.plot_complexity_vs_var(complexity_measure=cm, xvar='avg_temp')