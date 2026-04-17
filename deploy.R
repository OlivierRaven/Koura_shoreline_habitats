# deploy.R - runs automatically after quarto render
#system("powershell -Command \"Copy-Item -Recurse -Force '_manuscript\\*' 'docs\\'\"")
system("powershell -Command \"New-Item -ItemType Directory -Force 'docs\\analysis_files\\figure-html' | Out-Null\"")
system("powershell -Command \"Copy-Item -Force '.quarto\\_freeze\\analysis\\figure-html\\*' 'docs\\analysis_files\\figure-html\\'\"")
system("powershell -Command \"if (Test-Path 'site_libs') { Copy-Item -Recurse -Force 'site_libs' 'docs\\' }\"")
system("powershell -Command \"New-Item -Force 'docs\\.nojekyll' -ItemType File | Out-Null\"")
system("powershell -Command \"Copy-Item -Recurse -Force 'outputs' 'docs\\'\"")

# copy figures from both possible locations
system("powershell -Command \"New-Item -ItemType Directory -Force 'docs\\analysis_files\\figure-html' | Out-Null\"")
system("powershell -Command \"Copy-Item -Force 'analysis_files\\figure-html\\*' 'docs\\analysis_files\\figure-html\\' -ErrorAction SilentlyContinue\"")
system("powershell -Command \"Copy-Item -Force 'analysis_files\\figure-ipynb\\*' 'docs\\analysis_files\\figure-html\\' -ErrorAction SilentlyContinue\"")

cat("Deploy complete - docs folder updated\n")