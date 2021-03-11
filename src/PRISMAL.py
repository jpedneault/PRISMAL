#-------------------------Import-------------------------------------
import pandas as pd
from scipy import linalg
pd.options.mode.chained_assignment = None  # default='warn'

#----------------------------Class-----------------------------------
class Prismal:
    """ PRISMAL (PRoospective Impact Scenario Modelling of Aluminium) is a prospective tool to
        calculate impacts [GHG, EQ, HH] of the production of 1 kg of aluminium.
        3 main parameters changes in time:
            1) Change of electricity mix (based on the results of the SSP)
            2) Change of electrolysis technology (inert anode or prebaked)
            3) Improvments of energy intensity
    """

    # Initializer / Instance attributes
    def __init__(self, name, a, d_in, d_elec_ssp,d_elec_al, d_alu_prod, elec_int_pb, elec_int_in, electricity_tech, impacts, impact_list, regions, year,production):
        """Initialie all data needed for PRISMAL

         Args:
        -----
            name: [string] the name of the scenario
            a: [df] Technological matrix of aluminium production
            D_in: [df] Market share of inert anode technology
            d_elec_ssp: [df] Market share of electricity technology according SSP
            d_elec_al: [df] Market share of electricity technology for aluminium smelter
            D_alu_prod: [df] Marketshare of aluminum producer per region
            Elec_int_pb: [df] Electricity consumption of electrolysis (prebaked anode) process per year
            Elec_int_in: [df] Electricity consumption of electrolysis (inert anode) process per year
            Electricity_tech: [Serie] Name of all electriciy tech
            Impacts: [Serie] Name of all imapct categories
            Regions: [Serie] Name of all regions
            Year: [Serie] Years of calculation


        """

        # Initial atributes
        self.name = name
        self.a = a
        self.d_in = d_in
        self.d_elec_ssp = d_elec_ssp
        self.d_elec_al = d_elec_al
        self.d_alu_prod = d_alu_prod
        self.elec_int_pb = elec_int_pb
        self.elec_int_in = elec_int_in
        self.electricity_tech = electricity_tech
        self.impacts = impacts
        self.regions = regions
        self.year = year
        self.impact_list = impact_list
        self.production = production

        # To calculate
        self.l = pd.DataFrame()
        self.i_baux = pd.Series()
        self.i_oalu = pd.Series()
        self.i_smelting = pd.DataFrame()
        self.i_smelting_in = pd.DataFrame()
        self.i_smelting_pb = pd.DataFrame()
        self.i_anode = pd.DataFrame()
        self.i_anode_in = pd.DataFrame()
        self.i_anode_pb = pd.DataFrame()
        self.i_elec_smelting = pd.DataFrame()
        self.i_elec_smelting_in = pd.DataFrame()
        self.i_elec_smelting_pb = pd.DataFrame()
        self.i_cast = pd.DataFrame()
        self.i_elec_kWh = pd.DataFrame()
        self.i_total = pd.DataFrame()
        self.i_improvments = pd.DataFrame()

    @classmethod
    def from_template(cls, template_name, scenario_name, mitigation, prod_geo):
        """ Import all data needed for the calculation from the Excel sheet.
        First, we read all the excel file and sheet needed
        Second, we select the SSP according the scenario name
        Thrid, we create a scenrio object with all the data

        Args:
        -----
            template_name: [string] the name of the excel file to read (PRISMAL_data.xlsx)
            scenario_name: [string] the name of the PRISMAL scenario
            mitigation = [string] the name of the mitigation scenario selected

        Returns:
        --------
            scenario_object<>
        """

        # read excel
        a = pd.read_excel(template_name, index_col=[0], sheet_name='a')
        d_in = pd.read_excel(template_name, index_col=[0, 1], sheet_name='d_in')
        d_elec_ssp = pd.read_excel(template_name, index_col=[0, 1, 2, 3], sheet_name='d_elec_ssp')
        d_elec_al = pd.read_excel(template_name, index_col=[0, 1, 2, 3], sheet_name='d_elec_al')
        d_alu_prod = pd.read_excel(template_name, index_col=[0, 1], sheet_name=prod_geo)
        elec_int_in = pd.read_excel(template_name, index_col=[0, 1], sheet_name='elec_int_in')
        elec_int_pb = pd.read_excel(template_name, index_col=[0, 1], sheet_name='elec_int_pb')
        impacts = pd.read_excel(template_name, index_col=[0], sheet_name='impacts')
        production = pd.read_excel(template_name, index_col=[0], sheet_name='prod')

        #lists used to build index
        impact_list = list(impacts.index)
        year = list(d_elec_ssp.columns)
        electricity_tech = list(d_elec_ssp.index.get_level_values(3).unique())
        regions = list(d_elec_ssp.index.get_level_values(2).unique())
        regions.append('GLO')

        # Selection of scenario
        d_in = d_in.loc[scenario_name]
        d_in.columns = year
        d_elec_ssp = d_elec_ssp.loc[scenario_name].loc[mitigation]
        d_elec_al = d_elec_al.loc[scenario_name].loc[mitigation]
        d_alu_prod = d_alu_prod.loc[scenario_name]
        elec_int_in = elec_int_in.loc[scenario_name]
        elec_int_pb = elec_int_pb.loc[scenario_name]
        production = production.loc[scenario_name]

        # create object
        scenario_object = cls(scenario_name, a, d_in,  d_elec_ssp, d_elec_al, d_alu_prod, elec_int_in, elec_int_pb,
                              electricity_tech, impacts, impact_list, regions, year, production)
        return scenario_object

    def solve(self):
        """Solve all methods

        """
        self._calc_d_elec_al()
        self._calc_l()
        self._calc_i_baux()
        self._calc_i_oalu()
        self._calc_i_elec_kWh()
        self._calc_i_elec_smelting_pb()
        self._calc_i_elec_smelting_in()
        self._calc_i_elec_smelting()
        self._calc_i_smelting_in()
        self._calc_i_smelting_pb()
        self._calc_i_smelting()
        self._calc_i_anode_pb()
        self._calc_i_anode_in()
        self._calc_i_anode()
        self._calc_i_cast()
        self._calc_i_al_prod()
        self._calc_i_improvments()

    def solve_sens_elec(self):
        """Solve all methods

        """
        self._calc_d_elec_al()
        self._calc_l()
        self._calc_i_baux()
        self._calc_i_oalu()
        self._calc_i_elec_smelting_pb()
        self._calc_i_elec_smelting_in()
        self._calc_i_elec_smelting()
        self._calc_i_smelting_in()
        self._calc_i_smelting_pb()
        self._calc_i_smelting()
        self._calc_i_anode_pb()
        self._calc_i_anode_in()
        self._calc_i_anode()
        self._calc_i_cast()
        self._calc_i_al_prod()
        self._calc_i_improvments()

    def _calc_d_elec_al(self):
        """ Calculate future mix electric of aluminium smelter based on SSP electricity mix trends

        Args:
        -----

        Returns:
        --------
            self.d_elec_al: [df]
        """


        for r in self.regions[0:5]:



            d_plus = self.d_elec_ssp.loc[r].diff(axis=1).fillna(0)
            d_plus[d_plus < 0] = 0
            d_minus = self.d_elec_ssp.loc[r].diff(axis=1).fillna(0)
            d_minus[d_minus > 0] = 0
            delta_minus = self.d_elec_ssp.loc[r].pct_change(axis=1).fillna(0)
            delta_minus[delta_minus > 0] = 0
            delta_plus = d_plus / d_plus.sum(axis=0)
            delta_plus = delta_plus.fillna(0)

            d_elec_al_minus = pd.DataFrame(index=self.d_elec_al.index, columns=self.d_elec_al.columns)
            d_elec_al_plus = pd.DataFrame(index=self.d_elec_al.index, columns=self.d_elec_al.columns)

            for y in self.year:
                if y == 2020:
                    y_0 = y
                elif y == 2010:
                    y_0 = y
                elif y == 2005:
                    y_0 = y
                else:
                    d_elec_al_minus[y].loc[r] = self.d_elec_al[y_0].loc[r].values * delta_minus[y].values
                    d_elec_al_plus[y].loc[r] = -d_elec_al_minus[y].loc[r].sum() * delta_plus[y].values
                    self.d_elec_al[y].loc[r] = self.d_elec_al[y_0].loc[r].values + d_elec_al_minus[y].loc[r].values + d_elec_al_plus[y].loc[r].values
                    y_0 = y

    def _calc_l(self):
        """ Calculate inverse of the A matrix

        Args:
        -----

        Returns:
        --------
            self.l: [df]
        """

        self.l = pd.DataFrame(linalg.inv(self.a), index=self.a.index, columns=self.a.columns)

    def _calc_i_baux(self):
        """ Calculate the impact of bauxite needed for producing 1 kg of aluminium

        Args:
        -----

        Returns:
        --------
            self.i_baux: [df] impact per categories for every years
        """
        self.i_baux = pd.DataFrame(self.l.loc['bauxite', 'alu_ing'] * self.impacts.loc[:, 'Bauxite'])
        self.i_baux = pd.concat( [pd.concat([self.i_baux] * len(self.year), axis=1)] * len(self.regions), axis=0)

        stage = ['Bauxite']
        muix_baux = pd.MultiIndex.from_product([stage, self.regions, self.impact_list])
        self.i_baux.index = muix_baux
        self.i_baux.index.names = ['Stage', 'Regions', 'Impacts']
        self.i_baux.columns = self.year


    def _calc_i_oalu(self):
        """ Calculate the impact of alumina (o_alu) needed for producing 1 kg of aluminium

        Args:
        -----

        Returns:
        --------
            self.i_al_g: [df] impact per categories
        """

        self.i_oalu = pd.DataFrame(self.l.loc['o_alu', 'alu_ing'] * self.impacts.loc[:, 'Alumina'])
        self.i_oalu = pd.concat( [pd.concat([self.i_oalu] * len(self.year), axis=1)] * len(self.regions), axis=0)

        stage = ['Alumina']
        muix_oalu = pd.MultiIndex.from_product([stage, self.regions, self.impact_list])
        self.i_oalu.index = muix_oalu
        self.i_oalu.index.names = ['Stage', 'Regions', 'Impacts']
        self.i_oalu.columns = self.year

    def _calc_i_elec_kWh(self):
        """ Calculate the impact of 1 kWh of electricity per region and per year
        The method read the impact for every type of technology, and multiply them to the market share of
        every technology for every year (based on SSP results). To do the matrix multiplication, a diagonalisation
        of the impac is needed first

        Args:
        -----

        Returns:
        --------
            self.i_elec_kWh: [df] impact of electricity for every year
        """

        imp_elec = self.impacts.loc[:, self.electricity_tech]

        imp_elec_diag = pd.DataFrame(linalg.block_diag(*[imp_elec] * len(self.regions[0:5])),
                                     index=pd.MultiIndex.from_product([self.regions[0:5], imp_elec.index]),
                                     columns=self.d_elec_al.index)

        self.i_elec_kWh = imp_elec_diag @ self.d_elec_al

    def _calc_i_elec_smelting(self):
        """ Calculate the impact of electricity for the smelting process

        Args:
        -----

        Returns:
        --------
            self.i_elec_smelting: [df] impact of electricity for smelting
        """
        stage = ['Electricity for smelting']

        # Add electricity impact from Inert anode and prebaked anode
        self.i_elec_smelting = self.i_elec_smelting_in.mul(self.d_in, axis=0, level=0) + self.i_elec_smelting_pb.mul(1 - self.d_in, axis=0, level=0)

        # Calculate global impact with weigthd average
        self.d_alu_prod.columns = self.year
        i_elec_smelting_g = self.i_elec_smelting.mul(self.d_alu_prod, axis=0, level=0, fill_value=0).groupby(level=1).sum()
        i_elec_smelting_g.columns = self.year
        i_elec_smelting_g = pd.concat([i_elec_smelting_g], keys=['GLO'])

        # Concat region result and global result
        self.i_elec_smelting = pd.concat([self.i_elec_smelting, i_elec_smelting_g], axis=0)
        self.i_elec_smelting = pd.concat([self.i_elec_smelting], keys=stage)

        self.i_elec_smelting.index.names = ['Stage', 'Regions', 'Impacts']

    def _calc_i_elec_smelting_pb(self):
        """ Calculate the impact of electricity for the smelting process of prebaked anode

        Args:
        -----

        Returns:
        --------
            self.i_elec_smelting_pb: [df] impact of electricity for smelting
        """
        self.i_elec_smelting_pb = self.i_elec_kWh.mul(self.elec_int_pb, level=0)
        self.i_elec_smelting_pb.columns = self.year

    def _calc_i_elec_smelting_in(self):
        """ Calculate the impact of electricity for the smelting process of prebaked anode

        Args:
        -----

        Returns:
        --------
            self.i_elec_smelting_pb: [df] impact of electricity for smelting
        """

        self.i_elec_smelting_in = self.i_elec_kWh.mul(self.elec_int_in, level=0)
        self.i_elec_smelting_in.columns = self.year

    def _calc_i_smelting(self):
        """ Calculate the impact of smelting process without anode production and electricity

        Args:
        -----

        Returns:
        --------
            self.i_smelting: [df] impact of smelting
        """
        stage = ['Smelting']

        #Create matrix for every year based on impacts of smelting
        i_smelting_pb_s = pd.concat( [pd.concat([self.i_smelting_pb] * len(self.year), axis=1)] * len(self.regions[0:5]), axis=0)
        i_smelting_in_s = pd.concat( [pd.concat([self.i_smelting_in] * len(self.year), axis=1)] * len(self.regions[0:5]), axis=0)

        muix = pd.MultiIndex.from_product([self.regions[0:5], self.impact_list])
        i_smelting_pb_s.index = muix
        i_smelting_pb_s.columns = self.year
        i_smelting_in_s.index = muix
        i_smelting_in_s.columns = self.year

        #Impact of smelting for in and pb according market part
        i_smelting_in_d = i_smelting_in_s.mul(self.d_in, axis=0, level=0)
        i_smelting_pb_d = i_smelting_pb_s.mul(1 - self.d_in, axis=0, level=0)
        self.i_smelting = i_smelting_in_d + i_smelting_pb_d
        self.i_smelting.columns = self.year

        # Calculate global impact with weigthd average

        i_smelting_g = self.i_smelting.mul(self.d_alu_prod, axis=0, level=0, fill_value=0).groupby(level=1).sum()
        i_smelting_g = pd.concat([i_smelting_g], keys=['GLO'])
        i_smelting_g.columns = self.year

        # Concat region result and global result

        self.i_smelting = pd.concat([self.i_smelting, i_smelting_g], axis=0)
        self.i_smelting = pd.concat([self.i_smelting], keys=stage)
        self.i_smelting.index.names = ['Stage', 'Regions', 'Impacts']

    def _calc_i_smelting_pb(self):
        """ Calculate the impact of smelting process of prebake anode process without anode production and electricity

        Args:
        -----

        Returns:
        --------
            self.i_smelting: [df] impact of smelting
        """
        self.i_smelting_pb = pd.DataFrame(self.l.loc['alu_liq', 'alu_ing'] * self.impacts.loc[:, 'Smelting_pb'])

    def _calc_i_smelting_in(self):
        """ Calculate the impact of smelting process of inert anode process without anode production and electricity

        Args:
        -----

        Returns:
        --------
            self.i_smelting: [df] impact of smelting
        """
        self.i_smelting_in = pd.DataFrame(self.l.loc['alu_liq', 'alu_ing'] * self.impacts.loc[:, 'Smelting_in'])

    def _calc_i_anode(self):
            """ Calculate the impact of anode production by agregating prebaked and inert anode results

            Args:
            -----

            Returns:
            --------
                self.i_anode: [df] impact of smelting
            """

            stage = ['Anode']
            muix = pd.MultiIndex.from_product([self.regions[0:5], self.impact_list])

            # Create matrix for every year based on impacts of smelting

            i_anode_pb_s = pd.concat([pd.concat([self.i_anode_pb] * len(self.year), axis=1)] * len(self.regions[0:5]), axis=0)
            i_anode_in_s = pd.concat( [pd.concat([self.i_anode_in] * len(self.year), axis=1)] * len(self.regions[0:5]), axis=0)

            i_anode_pb_s.index = muix
            i_anode_pb_s.columns = self.year
            i_anode_in_s.index = muix
            i_anode_in_s.columns = self.year

            # Impact of anode production for in and pb according market part

            self.i_anode = i_anode_in_s.mul(self.d_in, axis=0, level=0) + i_anode_pb_s.mul(1 - self.d_in, axis=0, level=0)

            # Calculate global impact with weighted average

            i_anode_g = self.i_anode.mul(self.d_alu_prod, axis=0, level=0, fill_value=0).groupby(level=1).sum()
            i_anode_g = pd.concat([i_anode_g], keys=['GLO'])

            # Concat region result and global result

            self.i_anode = pd.concat([self.i_anode, i_anode_g], axis=0)
            self.i_anode = pd.concat([self.i_anode], keys=stage)
            self.i_anode.index.names = ['Stage', 'Regions', 'Impacts']

    def _calc_i_anode_pb(self):
            """ Calculate the impact of anode prebaked anode production

            Args:
            -----

            Returns:
            --------
                self.i_elec_smelting_pb: [df] impact of smelting """

            self.i_anode_pb = pd.DataFrame(self.l.loc['anode_pb', 'alu_ing'] * self.impacts.loc[:, 'Anode_pb'])

    def _calc_i_anode_in(self):
            """ Calculate the impact of inert anode production

            Args:
            -----

            Returns:
            --------
                self.i_elec_smelting_in: [df] impact of smelting """

            self.i_anode_in = pd.DataFrame(self.l.loc['anode_in', 'alu_ing'] * self.impacts.loc[:, 'Anode_in'])

    def _calc_i_cast(self):
        """ Calculate the impact of casting aluminium in ingot from liquid aluminium with changing electricity mix

        Args:
        -----

        Returns:
        --------
            self.i_cast: [df] impact of casting for every year
        """

        imp_casting = pd.DataFrame(self.l.loc['alu_ing', 'alu_ing'] * self.impacts.loc[:, 'Casting'])

        # Impact electricity for casting
        imp_cast_elec = self.i_elec_kWh.mul(self.l.loc['elec_cast', 'alu_ing'])
        self.i_cast = imp_cast_elec.add(imp_casting.squeeze(), axis=0, level=1)

        stage  = ['Casting']
        muix_cast = pd.MultiIndex.from_product([self.regions[0:5], self.impact_list])
        self.i_cast.index = muix_cast
        self.i_cast.columns = self.year

        # Calculate global impact with weighted average

        i_cast_g = self.i_cast.mul(self.d_alu_prod, axis=0, level=0, fill_value=0).groupby(level=1).sum()
        i_cast_g.columns = self.year
        i_cast_g = pd.concat([i_cast_g], keys=['GLO'])

        # Concat region result and global result

        self.i_cast = pd.concat([self.i_cast, i_cast_g], axis=0)
        self.i_cast = pd.concat([self.i_cast], keys=stage)
        self.i_cast.index.names = ['Stage', 'Regions', 'Impacts']

    def _calc_i_al_prod(self):
        """ Concatanation of every life cycle stages

        Args:
        -----

        Returns:
        --------
            self.i_cast: [df] impact of casting for every year
        """

        self.i_al_prod = pd.concat([self.i_baux, self.i_oalu, self.i_anode, self.i_elec_smelting, self.i_smelting, self.i_cast], axis=0)

    def _calc_i_improvments(self):
        """ Calculation fo the contribution of each improvments
        1) Changes in electricity mix
        2) Energy intensity improvments
        3) Inert anode deployment

        Args:
        -----

        Returns:
        --------
            self.i_improvments: [df] reduction of impact due to 3 parameters
        """


        # Energy intensity improvments

        first_year = self.year[2]

        elec_int = self.elec_int_in.mul(self.d_in, axis=0, level=0) + self.elec_int_pb.mul(1 - self.d_in, axis=0, level=0)
        elec_int_g = self.elec_int_in.mul(self.d_in, axis=0, level=0).mul(self.d_alu_prod, axis=0, level=0,
                                                                                  fill_value=0).sum(
            axis=0) + self.elec_int_pb.mul(1 - self.d_in, axis=0, level=0).mul(self.d_alu_prod, axis=0,
                                                                                       level=0, fill_value=0).sum(axis=0)

        delta_int = -elec_int.sub(elec_int.iloc(axis=1)[2], axis=0)
        delta_int_g = -elec_int_g.sub(elec_int_g.iloc[2])
        delta_int_g = pd.DataFrame(delta_int_g).T
        delta_int_g.index = ['GLO']

        i_elec = self.i_elec_kWh.loc(axis=1)[first_year]
        i_elec_g = i_elec.mul(self.d_alu_prod[first_year], axis=0, level=0).groupby(level=1).sum()
        i_elec_g = pd.concat([i_elec_g], keys=['GLO'])

        imp_ei = delta_int.mul(i_elec, axis=0, level=0)
        imp_ei_g = delta_int_g.mul(i_elec_g, axis=0, level=0)

        imp_ei = pd.concat([imp_ei, imp_ei_g], axis=0)
        imp_ei = pd.concat([imp_ei], keys=['Energy int. imp.'])
        imp_ei.index.names = ['Improvments', 'Regions', 'Impacts']

        # Changes in Electricity mix
        delta_i_elec = -self.i_elec_kWh.sub(self.i_elec_kWh[first_year], axis=0)
        i_elec_kwh_g = self.i_elec_kWh.mul(self.d_alu_prod, axis=0, level=0, fill_value=0).groupby(
            level=1).sum()
        delta_i_elec_g = -i_elec_kwh_g.sub(i_elec_kwh_g[first_year], axis=0)
        delta_i_elec_g = pd.concat([delta_i_elec_g], keys=['GLO'])

        imp_em = delta_i_elec.mul(elec_int, axis=0, level=0)
        elec_int_g = elec_int.mul(self.d_alu_prod, axis=0, level=0, fill_value=0).sum(axis=0)
        imp_em_g = delta_i_elec_g.mul((elec_int_g), axis=1, level=0)

        imp_em = pd.concat([imp_em, imp_em_g], axis=0)
        imp_em = pd.concat([imp_em], keys=['Changes in elec. mix'])
        imp_em.index.names = ['Improvments', 'Regions', 'Impacts']

        # Inert anode deployment

        muix = pd.MultiIndex.from_product([self.regions[0:5], self.impact_list])

        i_anode_pb = pd.concat([pd.concat([self.i_anode_pb] * len(self.year), axis=1)] * len(self.regions[0:5]), axis=0)
        i_anode_pb.index = muix
        i_anode_pb.columns = self.year

        i_smelting_pb = pd.concat([pd.concat([self.i_smelting_pb] * len(self.year), axis=1)] * len(self.regions[0:5]), axis=0)
        i_smelting_pb.index = muix
        i_smelting_pb.columns = self.year

        i_anode_in = pd.concat([pd.concat([self.i_anode_in] * len(self.year), axis=1)] * len(self.regions[0:5]), axis=0)
        i_anode_in.index = muix
        i_anode_in.columns = self.year

        i_smelting_in = pd.concat([pd.concat([self.i_smelting_in] * len(self.year), axis=1)] * len(self.regions[0:5]), axis=0)
        i_smelting_in.index = muix
        i_smelting_in.columns = self.year

        delta_smelting_tech = (i_anode_pb + i_smelting_pb - i_anode_in - i_smelting_in)
        imp_in = delta_smelting_tech.mul(self.d_in, axis=0, level=0)

        imp_in_g = imp_in.mul(self.d_alu_prod, axis=0, level=0, fill_value=0).groupby(level=1).sum()
        imp_in_g = pd.concat([imp_in_g], keys=['GLO'])

        imp_in = pd.concat([imp_in, imp_in_g], axis=0)
        imp_in = pd.concat([imp_in], keys=['Inert anode deployment'])
        imp_in.index.names = ['Improvments', 'Regions', 'Impacts']

        #Concatenate 3 improvments
        self.i_improvments = pd.concat([imp_ei, imp_em, imp_in], axis=0)

