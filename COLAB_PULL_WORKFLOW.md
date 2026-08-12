# Colab Pull Workflow

Code is stored in GitHub. Data and analysis outputs stay in Google Drive.
After any fresh Colab runtime, mount Drive and clone or update the repository:

```python
from google.colab import drive
drive.mount('/content/drive', force_remount=False)

import subprocess
from pathlib import Path

REPO_URL = 'https://github.com/anairr18/gastric_tme_project.git'
REPO = Path('/content/gastric_tme_project')
if REPO.exists():
    subprocess.run(['git', '-C', str(REPO), 'pull', '--ff-only'], check=True)
else:
    subprocess.run(['git', 'clone', REPO_URL, str(REPO)], check=True)
```

To run native per-cohort state discovery after the Korean ID repair has already
completed, use:

```python
import subprocess, sys
subprocess.run([
    sys.executable,
    '/content/gastric_tme_project/colab/run_native_per_cohort_state_discovery_colab.py',
], check=True)
```

This runner uses Drive inputs at `MyDrive/data`, writes results below
`MyDrive/gastric_tme_project`, and produces a review workbook. It never commits
or uploads any input data to GitHub.
