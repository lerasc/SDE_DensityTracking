from setuptools import setup

setup(
    name            =  'SDE_DensityTracking',
    version         =  '1.0.0',
    packages        = ['SDE_DensityTracking'],
    url             =  '',
    license         =  '',
    author          ='Sandro C. Lera',
    author_email    ='sandrolera@gmail.com',
    description     ='numerical solution of stochastic differential equation',
    python_requires ='>3.5.2',
    install_requires=[
                        "numpy>=1.20.3",
                        "pandas>=1.3.4",
                        "seaborn>=0.11.2",
                        "scipy>=1.7.2",
                     ]
)
