
CPPFLAGS = -I./
CXXFLAGS = -pipe -O3 -fshow-column -ffast-math -DEIGEN_DONT_VECTORIZE -fPIC \
	-DEIGEN_DONT_ALIGN -DNDEBUG -g -ffunction-sections -fdata-sections -DTIME_OPS
CXXFLAGS += -Wall -Wextra
LDFLAGS = -pg -g -L. -Wl,--no-undefined
LIBS = -lboost_thread

#ifneq ($(WINDIR),)
# On windows
#EXEXT = .exe
#CXX = mingw32-g++ -std=gnu++0x
#CC = mingw32-gcc -std=c89
#LD = $(CC)
#CXXFLAGS += -mincoming-stack-boundary=2 -march=i686 -msse2
#CPPFLAGS += -I./win32
#CFLAGS = -O1 -g
#TARGETDIR = win32
#else
#assumed to be on Linux
EXEEXT = 
CXX = g++-4.4 -std=gnu++0x
CC = gcc-4.4
LD = $(CC)
LDXX = $(CXX)
CPPFLAGS += -I./posix -I/home/jonathan/src/eigen -I/home/jonathan/src/boost_1_43_0
TARGETDIR = posix
#endif

LDXX = $(CXX)
RM = rm
AR = ar

all:  test_ins_qkf$(EXEXT) libeknav-ecef.so libeknav-ned.so

tags : *.cpp
	ctags ./*.cpp ./*.hpp

clean:
	-rm *.o
	-rm $(TARGETDIR)/*.o
	-rm test_ins_qkf$(EXEXT)

INS_QKF_OBJS = ins_qkf_observe_gps_pvt.o \
	ins_qkf_predict_ned.o \
	ins_qkf_observe_gps_p.o \
	ins_qkf_observe_vector.o \
	ins_qkf_predict.o \
	diagnostics.o \
	basic_ins_qkf.o

INS_QKF_NED_OBJS = ins_qkf_observe_gps_p.o \
	ins_qkf_observe_vector.o \
	ins_qkf_predict_ned.o \
	basic_ins_qkf.o \
	$(TARGETDIR)/timer.o

INS_QKF_ECEF_OBJS = ins_qkf_observe_gps_pvt.o \
	ins_qkf_observe_gps_p.o \
	ins_qkf_observe_vector.o \
	ins_qkf_predict.o \
	basic_ins_qkf.o \
	$(TARGETDIR)/timer.o

PLATFORM_OBJS = $(TARGETDIR)/timer.o \
	$(TARGETDIR)/random_seed.o

%.o: %.cpp
	$(CXX) -c -o $@ $< -MMD $(CXXFLAGS) $(CPPFLAGS)

libeknav-ned.so: $(INS_QKF_NED_OBJS)
	$(LDXX) -shared -O3 -o $@ $^ $(LDFLAGS)

libeknav-ecef.so: $(INS_QKF_ECEF_OBJS)
	$(LDXX) -shared -O3 -o $@ $^ $(LDFLAGS)

libeknav-ecef.a: $(INS_QKF_ECEF_OBJS)
	$(AR) cru $@ $^

libeknav-ned.a: $(INS_QKF_NED_OBJS)
	$(AR) cru $@ $^

test_ins_qkf$(EXEXT): test_ins_qkf.o diagnostics.o $(PLATFORM_OBJS) libeknav-ecef.so
	$(LDXX) -o $@ $< diagnostics.o $(PLATFORM_OBJS) $(LDFLAGS) $(LIBS) -L. -leknav-ecef

-include *.d win32/*.d posix/*.d
