# inverse-participation-ratio-optimization-to-find-confined-states-in-quasicrystals
強束縛模型ハミルトニアンの対角化によって得られる準結晶の縮退状態をコンファインド状態へ構成しなおすコード。

## 使い方
### 1. coordinates.txt, hop.txt にそれぞれ準結晶の座標とホッピングの情報を入力する。
### 2. TB.f90 の N と allhop にそれぞれサイト数、ボンド数を代入し、make_Hamiltonian サブルーチンの edge_potential 部分を適当に決める。
### 3. TB.f90 を実行し、wave_function フォルダに波動関数のファイルを出力する。同時に生成された eigval ファイルを見ながら、縮退状態が適切に得られていることを確認したら、wave_function/AB{N} というようにフォルダを作成し、波動関数のファイルをそのフォルダに入れる。このとき、wf_plot.py を実行して波動関数の分布を見ることができる。
### 4. IPR.f90 ファイルの N, N_c にサイト数、縮退状態数を代入する。number を1に、is_to_update を .false. にして、以下の操作を開始する。
1. psis の初期値を適当に決め、IPR.f90 を実行し、confined_states/AB{N} に psis_{number}.txt を出力する。
2. 標準出力されたIPRの差分が十分小さく、それまでに得られたコンファインド状態と重複しないことを確認したら number を 1 つ上げ、is_to_update を .false. にして 1. を再び行う。そうでなければ、is_to_update を .true. にして 1. を行う。
3. number が N_c より大きくなったら終了する。
### 5. confined_states_plot.py を実行し、コンファインド状態の波動関数強度のプロットを得る。
