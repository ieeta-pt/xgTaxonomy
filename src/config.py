from pathlib import Path

base_path = Path(__file__).parent
summariesPath = (base_path / "../aux/summaries").resolve()
dstPathOfOriginalSequences = (base_path / "../aux/original_sequences").resolve()
dstPathOfOriginalAminoacids = (base_path / "../aux/original_protein").resolve()
locationOfDatabases = (base_path / "../aux/references").resolve()

dstPathOfDatabaseSequences = (base_path / "../aux/database_sequences").resolve()
dstPathOfDatabaseAminoacids = (base_path / "../aux/database_protein").resolve()
orfPath = (base_path / "../ORFs/orfm/").resolve()
genomeFeaturesFilePath = (base_path / "../results/genome_features.csv").resolve()
proteomeFeaturesFilePath = (base_path / "../results/proteome_features.csv").resolve()
tmpPrt=(base_path / "../aux/prt").resolve()
numIterations = 10