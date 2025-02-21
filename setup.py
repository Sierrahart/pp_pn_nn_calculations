from setuptools import setup

setup(
    name='PP_conformers',
    version='0',
    packages=['pp_conformers'],
    url='',
    license='',
    author='Jamie A. Cadge and Sierra Hart',
    author_email='j.cadge@icloud.com',
    description='',
    python_requires=">=3.8",
    install_requires=["mdanalysis==2.2.0", "morfeus-ml==0.7.2",
                      "numpy==1.24.4", "pandas==2.0.3",
                       "yaml"]
)
