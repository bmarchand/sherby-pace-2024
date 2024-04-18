import sys
import os

TEST_VS_KOBA = 0
TEST_EXCLUDED = 1

mode = TEST_VS_KOBA

exclude = [31, 84] + list(range(38, 55)) + list(range(58, 70)) + list(range(73, 83)) + list(range(92, 100))


#this erases all contents of file
open("test.txt", 'w').close()

if mode == TEST_EXCLUDED:
    for i in exclude:
        cmd = "timeout 30m cargo run --release --bin=sherby-pace-2024 -- exact-public-instances/" + str(i) + ".gr --dfas >> test.txt"
        print(cmd)
        os.system(cmd)

elif mode == TEST_VS_KOBA:
    



    for i in range(1, 101):
        if not i in exclude:    #suboptimal O(n)
            
            cmd = f'echo -e "\\n{i}.gr" >> test.txt'   #lol at this way of outputting into test file
            os.system(cmd)
            
            cmd = "cargo run --release --bin=sherby-pace-2024 -- exact-public-instances/" + str(i) + ".gr"
            print(cmd)
            os.system(cmd)
            
            
            cmd = "pace2024verifier exact-public-instances/" + str(i) + ".gr exact-public-instances/" + str(i) + ".sol -c >> test.txt"
            print(cmd)
            os.system(cmd)
            
            #just to be sure we don't reevaluate the same file
            cmd = "rm exact-public-instances/" + str(i) + ".sol"
            print(cmd)
            os.system(cmd)
            
            cmd = "cargo run --release --bin=sherby-pace-2024 -- exact-public-instances/" + str(i) + ".gr --dfas"
            print(cmd)
            os.system(cmd)
            
            
            cmd = "pace2024verifier exact-public-instances/" + str(i) + ".gr exact-public-instances/" + str(i) + ".sol -c >> test.txt"
            print(cmd)
            os.system(cmd)
            