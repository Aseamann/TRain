from setuptools import setup, find_packages

# To perform build type "python -m pip install -e ."

setup(
    name="TRain",
    version="0.1.0",
    description="TRain is a tool that...",
    author="Austin Seamann",
    author_email="aseamann@unomaha.edu",
    url="https://github.com/Aseamann/TRain",
    install_requires=[
        "pandas",
        "openpyxl",
        "biopython == 1.79",
        "scipy",
        "scikit-learn",
        "matplotlib",
        "seaborn"
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