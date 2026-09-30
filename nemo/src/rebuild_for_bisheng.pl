INUNITS=( numnam_ref numnam_cfg numnat_ref numnat_cfg numtrc_ref numtrc_cfg numnam_ice_ref numnam_ice_cfg  numnatp_cfg numnatp_ref )

# 构建需要修改的文件列表
# 查找与 INUNITS 匹配的文件
for iunit in "${INUNITS[@]}"; do
  grep -l "$iunit" `find ./ -name '*.[fFh]90'` >> temp_list.txt
done

# 删除重复项并生成最终文件列表
allfiles=$(sort -u temp_list.txt)

# 删除临时文件
rm temp_list.txt

# 对每个找到的文件应用 Perl 替换
for f in $allfiles; do
  echo "Working on $f"

  # 对于 INUNITS 中的每个单位
  for n in `seq 0 1 $(( ${#INUNITS[*]} - 1 ))`; do
    numnam=${INUNITS[$n]}

    # 根据条件执行替换
    perl -pi -e '
    $pattern = qr/(.*\s*READ\s*\()(\s*)('${INUNITS[$n]}')(\s*,\s*)([a-z0-9_]*)(.*)/i;

    if ($_ =~ $pattern) {
      # 如果第五组是 "namtrc"，则在 INDEX 中增加空格
      if ($5 eq "namtrc") {
        s@$pattern@$1$2$3(INDEX($3,"$5 ")-1:)$4$5$6@;
      } else {
        # 否则，保持原逻辑
        s@$pattern@$1$2$3(INDEX($3,"$5")-1:)$4$5$6@;
      }
    }
    ' "$f"
  done
done

