"""Shared project credits for repository guides and Zenodo metadata."""
import html

CREDITS = [
    dict(name='pycalphad', role='Thermodynamic calculations and phase equilibria',
         papers=[('R. Otis and Z.-K. Liu (2017), pycalphad: CALPHAD-based Computational Thermodynamics in Python', '10.5334/jors.140')],
         repository='https://github.com/pycalphad/pycalphad'),
    dict(name='ESPEI', role='Thermodynamic parameter fitting and uncertainty quantification',
         papers=[('B. Bocklund et al. (2019), ESPEI for efficient thermodynamic database development, modification, and uncertainty quantification: application to Cu–Mg', '10.1557/mrc.2019.59')],
         repository='https://github.com/PhasesResearchLab/ESPEI'),
    dict(name='PhaseForge', role='MLIP-informed CALPHAD construction in the examples',
         papers=[('S. Zhu, D. Sarıtürk, and R. Arróyave (2025), Machine learning potentials for alloys: a detailed workflow to predict phase diagrams and benchmark accuracy', '10.1038/s41524-025-01814-z')],
         repository='https://github.com/dogusariturk/PhaseForge', software_doi='10.5281/zenodo.15730912'),
    dict(name='PhaseForge+', role='Construction and experimental tuning of CALPHAD models in the examples; repository name PhaseForgePlus',
         papers=[('C. Kunselman, S. Zhu, D. Sarıtürk, and R. Arróyave (2025), Construction and Tuning of CALPHAD Models Using Machine-Learned Interatomic Potentials and Experimental Data: A Case Study of the Pt–W System', '10.1007/s11669-025-01222-2')],
         repository='https://github.com/dogusariturk/PhaseForgePlus', software_doi='10.5281/zenodo.17519087'),
    dict(name='Materials Project', role='DFT structures and energies; repository link is the API client, not the database itself',
         papers=[('A. Jain et al. (2013), Commentary: The Materials Project: A materials genome approach to accelerating materials innovation', '10.1063/1.4812323')],
         repository='https://github.com/materialsproject/api', website='https://materialsproject.org'),
    dict(name='OQMD', role='DFT structures and energies; qmpy is the supporting software',
         papers=[('J. E. Saal et al. (2013), Materials Design and Discovery with High-Throughput Density Functional Theory: The Open Quantum Materials Database (OQMD)', '10.1007/s11837-013-0755-4'),
                 ('S. Kirklin et al. (2015), The Open Quantum Materials Database (OQMD): assessing the accuracy of DFT formation energies', '10.1038/npjcompumats.2015.10')],
         repository='https://github.com/wolverton-research-group/qmpy', website='https://oqmd.org'),
]
INTRO = ('These credits describe the ORCHID project and its examples. Cite the papers and software used in your analysis, '
         'together with the ORCHID release and original source of each observation. The list does not imply that every '
         'dataset was generated with every package. Preserve the versions and commit hashes recorded by each example; '
         'a repository link is not a reproducibility pin. Materials Project and OQMD energies retain separate references.')


def markdown():
    lines = ['# Software and data citations', '', INTRO, '']
    for c in CREDITS:
        lines += ['## '+c['name'], '', c['role']+'.', '']
        lines += ['- '+title+'. [Paper](https://doi.org/'+doi+').' for title,doi in c['papers']]
        lines += ['- [Repository]('+c['repository']+').']
        if c.get('software_doi'):lines += ['- [Archived software](https://doi.org/'+c['software_doi']+'). Cite the version actually used.']
        if c.get('website'):lines += ['- [Data resource]('+c['website']+').']
        lines += ['']
    lines += ['The Cu–Zr gradient workflows additionally credit Courtney Kunselman’s ESPEI and littlemcmc branches, '
              'as pinned in the example requirements and provenance manifests: '
              '[analytic gradients](https://doi.org/10.1016/j.actamat.2025.121453), '
              '[gradient-informed HMC](https://doi.org/10.1016/j.actamat.2026.122412), '
              '[ESPEI fork](https://github.com/cjkunselman18/ESPEI), '
              '[littlemcmc fork](https://github.com/cjkunselman18/littlemcmc).', '']
    return '\n'.join(lines)


def html_block():
    h = '<section id="software-data-credits"><h2>Software and data citations</h2><p>'+html.escape(INTRO)+'</p><ul>'
    for c in CREDITS:
        h += '<li><b>'+html.escape(c['name'])+'</b>: '+html.escape(c['role'])+'. '
        h += '; '.join(html.escape(t)+' (<a href="https://doi.org/'+d+'">'+d+'</a>)' for t,d in c['papers'])
        h += '. <a href="'+c['repository']+'">Repository</a>'
        if c.get('software_doi'):h += '; <a href="https://doi.org/'+c['software_doi']+'">archived software</a>'
        if c.get('website'):h += '; <a href="'+c['website']+'">data resource</a>'
        h += '.</li>'
    return h+'</ul></section>'


def related_identifiers():
    out=[]
    for c in CREDITS:
        for _,doi in c['papers']:
            out.append(dict(identifier=doi,relation='references',scheme='doi'))
        out.append(dict(identifier=c['repository'],relation='references',scheme='url'))
        if c.get('software_doi'):
            out.append(dict(identifier=c['software_doi'],relation='references',scheme='doi',resource_type='software'))
    return out
