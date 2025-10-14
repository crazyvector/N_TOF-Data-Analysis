#!/bin/bash

# Locul de destinație (folderul unde vor fi copiate fișierele)
DEST="/eos/user/a/acristes/GitHub_Repositories/N_TOF-Data-Analysis/Runs"

# Creează folderul dacă nu există
mkdir -p "$DEST"

# Lista de run-uri
runs=(
119910 119911 119912 119913 119914 119915 119916 119917 119918 119919
119920 119921 119922 119923 119924 119926 119928 119929 119930 119933
119934 119935 119936 119937 119938 119939 119940 119941 119942 119946
119947 119948 119951 119953 119956 119957 119958 119959 119960 119972
119974 119975 119976 119977 119981
)

# Copiază toate fișierele din directorul oficial în spațiul tău EOS
for run in "${runs[@]}"; do
    echo "Copying run $run..."
    cp /eos/experiment/ntof/processing/official/done/run${run}.root "$DEST"/
done

echo "✅ All runs copied to $DEST"
