import numpy as np
import pandas as pd
from typing import List, Tuple
from networks import create_region_species_network
from complexity import compute_eci, compute_fitness, count
import plotly.graph_objects as go
from scipy.stats import linregress

class RegionData:
    """

    Wrapper for population density data

    """
    def __init__(self, geo_groups: List[str]):
        """

        :param geo_groups: the geographic region groups (e.g. Europe, China..) for which you want to load the data
        """
        self.geo_groups = geo_groups
        self.df_region_data = self._load_data(geo_groups=geo_groups)

    @staticmethod
    def _load_data(geo_groups: List[str]) -> pd.DataFrame:
        """

        This method loads population density data

        :param geo_groups: the geographic region groups (e.g. Europe, China..) for which you want to load the data
        :return: A dataframe with columns ['region', 'abbrev', 'density', 'macro_region'].

        """
        df = pd.DataFrame(columns=['region', 'forest_area', 'pop_density', 'macro_region', 'area', 'avg_temp', 'min_temp', 'max_temp'])
        for gg in geo_groups:
            name_csv = f'./region_data/{gg}Data.csv'
            df_mr = pd.read_csv(filepath_or_buffer=name_csv, header=0)
            df_mr['macro_region'] = gg
            df = pd.concat([df, df_mr], ignore_index=True)

        df.sort_values('region', inplace=True)
        return df

class EcosystemComplexityData:
    """

    Wrapper for economic complexity data of species and regions

    """
    def __init__(self, geo_groups: List[str]):
        """
        :param geo_groups: the geographic region groups (e.g. Europe, China..) for which you want to load the data
        """
        self.geo_groups = geo_groups
        self.df_regions, self.df_species = self._load_data(geo_groups=geo_groups)

    @staticmethod
    def _load_data(geo_groups: List[str]) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """

        Computes country and species complexity indices for all regions in the geographic groups and stores in a dataframe
        with columns ['regions' (or 'species'), 'cou_r', 'eci_r', 'fit_r'].

        :param geo_groups: the geographic region groups (e.g. Europe, China..) for which you want to load the data
        :return: A tuple of dataframes, one with complexity measures for regions and the other with complexity measures
        for species

        """
        tax_groups = []
        network, matrix, matrix_sorted, regions, regions_sorted, species, species_sorted, taxonomic_map, desc \
            = create_region_species_network(geo_groups, tax_groups)

        cou_s, cou_r = count(matrix=matrix)
        eci_s, eci_r = compute_eci(matrix=matrix, regions=regions, species=species)
        fit_s, fit_r = compute_fitness(matrix=matrix, regions=regions, species=species)

        data_regions = np.stack([regions, cou_r, eci_r, fit_r], axis=1)
        df_regions = pd.DataFrame(data=data_regions, columns=['region', 'cou_r', 'eci_r', 'fit_r'])
        df_regions.sort_values('region', inplace=True)

        data_species = np.stack([species, cou_s, eci_s, fit_s], axis=1)
        df_species = pd.DataFrame(data=data_species, columns=['species', 'cou_r', 'eci_r', 'fit_r'])
        df_species.sort_values('species', inplace=True)

        return df_regions, df_species


class Plotter:
    @staticmethod
    def plot(x: np.ndarray, y: np.ndarray, text: np.ndarray, xlabel: str, ylabel: str, xscale: str = 'log',
              yscale: str = 'log', marker_color: str = 'black', mode: str = 'markers', fit_line: bool = True,
              path_image: str = '', showscale: bool = False, showlegend: bool = True, title: str = ''):
        """

        Plotting function

        :param x: the x axis variable
        :param y: the y axis variable
        :param text: the text you want to appear when you hover over a point
        :param xlabel: the label of the x axis
        :param ylabel: the label of the y axis
        :param xscale: the type of the x axis, options are 'log' or 'linear'
        :param yscale: the type of the y axis, options are 'log' or 'linear'
        :param marker_color: the values with which to color the markers
        :param mode: variable to select scatter plot or line plot, options are 'markers' or 'lines'
        :param fit_line:  True if you want a line fit of your data and false otherwise
        :return: None

        """
        fig = go.Figure()
        trace = go.Scatter(x=x, y=y, text=text, mode=mode, name=ylabel,
                           marker=dict(color=marker_color, colorscale='Viridis', showscale=showscale))

        if fit_line:
            a, b, r = Plotter._fit(x=x, y=y, xscale=xscale, yscale=yscale)

            if xscale == 'log' and yscale == 'log':
                w = np.exp(b) * np.power(x, a)
            elif xscale == 'log' and yscale == 'linear':
                w = b + a * np.log(x)
            elif xscale == 'linear' and yscale == 'log':
                w = np.exp(a * x + b)
            elif xscale == 'linear' and yscale == 'linear':
                w = a * x + b
            else:
                w = np.zeros(len(x))

            trace_reg = go.Scatter(x=x, y=w, mode='lines',
                                   name=f'a={np.round(a, decimals=5)} b={np.round(b, decimals=5)}, R^2={np.round(r ** 2, decimals=5)}')
            fig.add_trace(trace_reg)

        fig.add_trace(trace)
        fig.update_xaxes(type=xscale)
        fig.update_yaxes(type=yscale)
        fig.update_layout(title=title,
                          xaxis_title=xlabel,
                          yaxis_title=ylabel,
                          showlegend=showlegend,
                          font=dict(size=18, color="Black"))

        if path_image == '':
            fig.show(renderer="browser")
        else:
            fig.write_image(path_image)

    @staticmethod
    def _fit(x: np.ndarray, y: np.ndarray, xscale: str, yscale: str) -> Tuple[float, float, float]:
        """

        Fits a line to your data

        :param x: the x axis variable
        :param y: the y axis variable
        :param xscale: the type of the x axis, options are 'log' or 'linear'
        :param yscale: the type of the y axis, options are 'log' or 'linear'
        :return: the slope, intercept and r of the fitted line (line is a*x + b)

        """
        if xscale == 'log' and yscale == 'log':
            v = np.log(x)
            w = np.log(y)
        elif xscale == 'log' and yscale == 'linear':
            v = np.log(x)
            w = y
        elif xscale == 'linear' and yscale == 'log':
            v = x
            w = np.log(y)
        elif xscale == 'linear' and yscale == 'linear':
            v = x
            w = y
        else:
            v = 0
            w = 0

        slope, intercept, r, p, se = linregress(v, w)
        return slope, intercept, r
