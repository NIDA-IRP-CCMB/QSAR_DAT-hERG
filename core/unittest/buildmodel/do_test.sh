#!/bin/bash
module load python


echo "running reg xgb"
python core/run_buildmodel.py -c regression -m xgb -i dataset_rDAT_inhibitor_Ki -b pubdata -o output_reg_xgb -t 0.00 -n 25 | tee output_reg_xgb.log
echo "check reg xgb"
diff output_reg_xgb.log reference/rDAT_inhibitor_reg_xgb_alldata.log

echo "running reg rf"
python core/run_buildmodel.py -c regression -m rf -i dataset_rDAT_inhibitor_Ki -b pubdata -o output_reg_rf -t 0.00 -n 25 | tee output_reg_rf.log
echo "check reg rf"
diff output_reg_rf.log reference/output_reg_rf.log

echo "running class xgb"
python core/run_buildmodel.py -c classification -m xgb -i dataset_rDAT_inhibitor_Ki -b pubdata_class -o output_class_xgb -t 0.00 -n 25 | tee output_class_xgb.log
echo "check class xgb"
diff output_class_xgb.log reference/rDAT_inhibitor_class_xgb_alldata.log

echo "running class rf"
python core/run_buildmodel.py -c classification -m rf -i dataset_rDAT_inhibitor_Ki -b pubdata_class -o output_class_rf -t 0.00 -n 25 | tee output_class_rf.log
echo "check class rf"
diff output_class_rf.log reference/rDAT_inhibitor_class_rf_alldata.log



