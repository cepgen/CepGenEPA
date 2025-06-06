import Config.Core as cepgen
from Config.PDG_cfi import PDG
from Config.timer_cfi import timer # enable timing framework
#from Integrators.miser_cfi import miser as integrator
#from Integrators.foam_cfi import foam as integrator


#integrator = cepgen.Module('bases')
from Config.logger_cfi import logger
#logger.enabledModules += ('Process.weight',)


process = cepgen.Module('epa',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.ElasticElastic,
        matrixElement = cepgen.Module('python',
            function = 'trapint_electron_Hamzeh.cs_electron_w_condition_Hamzeh',
        ),
        partonsFlux = cepgen.Module('grid',
            modelling = cepgen.Module('python',
                function = 'Integrated_elastic_tau_tau_cross_section_final_version.flux_el_yy_atW',
                beam1 = cepgen.Parameters(
                    energy = 50.,
                ),
                beam2 = cepgen.Parameters(
                    energy = 7000.,
                ),
            ),
            path = 'flux.grid',
            #generateGrid = True,  # force the grid (re-)computation
        ),
    ),
    inKinematics = cepgen.Parameters(
        pz = (50., 7000.),
        pdgIds = (11, 2212),
    ),
    outKinematics = cepgen.Parameters(
        invmass = (10.,),
    )
)

generator = cepgen.Parameters(
    numEvents = 10000
)
text = cepgen.Module('text',
    #variables = ['nev', 'm(4)', 'tgen'],
    histVariables={
        'm(4)': cepgen.Parameters(xbins=[float(bin) for bin in range(0, 250, 10)]),
        'm(ob1)': cepgen.Parameters(xrange=(0., 250.), nbins=10, log=True),
        'pt(7):pt(8)': cepgen.Parameters(xrange=(0., 250.), yrange=(0., 250.), log=True)
    }
)
output = cepgen.Sequence(text)
