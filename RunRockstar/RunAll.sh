./main.sh &
echo "Started Main Process"
sleep 10
./first.sh &
./second.sh &
./third.sh &
./fourth.sh
echo "Finished Rockstar Halo Finding"
./FinishRockstar.sh
