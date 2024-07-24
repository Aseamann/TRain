from setuptools import setup, find_packages

# To perform build type "python -m pip install -e ."
# Requires python==3.9

setup(
    name="TRain",
    version="1.0.0",
    description="TRain is a tool that allows users to prepare TCR docking data and perform docking.",
    author="Austin Seamann",
    author_email="aseamann@unomaha.edu",
    url="https://github.com/Aseamann/TRain",
    install_requires=[
        "pillow == 8.2.0"
        "pandas == 1.4.4",
        "numpy == 1.20.1",
        "openpyxl == 3.1.2",
        "biopython == 1.79",
        "scipy == 1.7.2",
        "scikit-learn == 1.0.1",
        "matplotlib == 3.7.4",
        "seaborn == 0.13.1",
        "bio == 1.5.0"
    ],
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'SeqConductor = tr01input:main',
            'ModelEngine = tr02model:main',
            'TurnTable = tr03pair:main',
            'TCRcoupler = tr04dock:main',
            'PrepCoupler = tr04dock:prep',
            'PostCoupler = tr04dock:post',
            'DataDepot = tr05analysis:main',
            'PDB_Tools_V3 = util:main'
        ]
    },
    include_package_data=True,
    data_files=[('config', ["data/config.ini"]), ('fasta', ["data/*.fasta"])],
    package_data={"tr03pair": ["*.pdb"]}
)
