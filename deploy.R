# deploy.R - runs automatically after quarto render

# Create the target directory for figures matching HTML paths
system("powershell -Command \"New-Item -ItemType Directory -Force 'docs\\index_files\\figure-html' | Out-Null\"")

# Copy figure files from Quarto freeze folder (if exists) to the target folder
system("powershell -Command \"Copy-Item -Force '.quarto\\_freeze\\analysis\\figure-html\\*' 'docs\\index_files\\figure-html\\' -ErrorAction SilentlyContinue\"")

# Copy figure files from other possible source folders to the target folder
system("powershell -Command \"Copy-Item -Force 'analysis_files\\figure-html\\*' 'docs\\index_files\\figure-html\\' -ErrorAction SilentlyContinue\"")
system("powershell -Command \"Copy-Item -Force 'analysis_files\\figure-ipynb\\*' 'docs\\index_files\\figure-html\\' -ErrorAction SilentlyContinue\"")

# Copy site libraries if they exist
system("powershell -Command \"if (Test-Path 'site_libs') { Copy-Item -Recurse -Force 'site_libs' 'docs\\' }\"")

# Copy outputs folder into docs
system("powershell -Command \"Copy-Item -Recurse -Force 'outputs' 'docs\\'\"")

# Create .nojekyll file to prevent GitHub Pages ignoring files/folders starting with _
system("powershell -Command \"New-Item -Force 'docs\\.nojekyll' -ItemType File | Out-Null\"")

cat("Deploy complete - docs folder updated\n")
