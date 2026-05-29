# deploy.R - post-render script: copies build artefacts into docs/ for GitHub Pages

# Ensure figure directories exist
system("powershell -Command \"New-Item -ItemType Directory -Force 'docs\\analysis_files\\figure-html' | Out-Null\"")
system("powershell -Command \"New-Item -ItemType Directory -Force 'docs\\index_files\\figure-html' | Out-Null\"")

# Copy figures from all possible source locations into docs/analysis_files/
system("powershell -Command \"Copy-Item -Force '.quarto\\_freeze\\analysis\\figure-html\\*' 'docs\\analysis_files\\figure-html\\' -ErrorAction SilentlyContinue\"")
system("powershell -Command \"Copy-Item -Force 'analysis_files\\figure-html\\*' 'docs\\analysis_files\\figure-html\\' -ErrorAction SilentlyContinue\"")
system("powershell -Command \"Copy-Item -Force 'analysis_files\\figure-ipynb\\*' 'docs\\analysis_files\\figure-html\\' -ErrorAction SilentlyContinue\"")

# Mirror into docs/index_files/ so index.html figure links resolve
system("powershell -Command \"Copy-Item -Force 'docs\\analysis_files\\figure-html\\*' 'docs\\index_files\\figure-html\\' -ErrorAction SilentlyContinue\"")

# Copy site libraries
system("powershell -Command \"if (Test-Path 'site_libs') { Copy-Item -Recurse -Force 'site_libs' 'docs\\' }\"")

# Copy outputs folder
system("powershell -Command \"Copy-Item -Recurse -Force 'outputs' 'docs\\'\"")

# Prevent GitHub Pages from ignoring underscore-prefixed folders
system("powershell -Command \"New-Item -Force 'docs\\.nojekyll' -ItemType File | Out-Null\"")

cat("Deploy complete - docs folder updated\n")
