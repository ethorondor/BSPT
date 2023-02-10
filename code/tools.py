class Computer:
    def __init__(self, cpu, ram):
        self.cpu = cpu
        self.ram = ram
        
    def config(self):
        print('config is', self.cpu, self.ram)

def hello_world(name):
    print("hello {0}".format(name))