#include <iostream>
#include <chrono>
#include <vector>
#include <cmath>

class timer
{
    class measure_point
    {
    public:
        std::string text;
        std::chrono::high_resolution_clock::time_point start;
        measure_point(std::string text) : text(text),start(std::chrono::high_resolution_clock::now()) {}
    };

    std::string text;
    std::vector<measure_point> points;
    bool output;            ///< if the timer should print something on destruction
public:
    std::chrono::high_resolution_clock::time_point start;
    /// \brief The text is the general name for this timer. Setting output to false disables printing anything on destruction.
    timer(std::string text="",bool output=true) : text(text),output(output) {start=std::chrono::high_resolution_clock::now();}
    /// \brief The text is the general name for this timer. Setting output to false disables printing anything on destruction.
    timer(const char* text,   bool output=true) : text(text),output(output) {start=std::chrono::high_resolution_clock::now();}
    timer(const timer&)=default;
    timer(timer&&)=default;
    timer& operator=(const timer&)=default;
    timer& operator=(timer&&)=default;

    ~timer()
    {
        if(!output)
            return;
        auto start=this->start;
        double diff=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-start).count()/1000000.0;
        if(text.size())
            std::cout<<diff<<" \t<- "<<text<<std::endl;
        else
            std::cout<<"timer time: "<<diff<<std::endl;

        if(points.empty())
            return;

        for(auto p:points)
        {
            auto diff=std::chrono::duration_cast<std::chrono::microseconds>(p.start-start).count()/1000000.0;
            std::cout<<"  "<<diff<<" \t<- "<<p.text<<std::endl;
            start=p.start;
        }
        diff=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-start).count()/1000000.0;
        std::cout<<"  "<<diff<<" to end"<<std::endl;
    }

    /// \brief Returns the time passed since the timer was started or reset in seconds.
    double until_now() const
    {
        return std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-start).count()/1000000.0;
    }

    void add(const std::string& name="")
    {
        if(!name.size())
            points.emplace_back(std::to_string(points.size()));
        else
            points.emplace_back(std::to_string(points.size())+" "+name);
    }

    /// \brief Resets the starting time to the current time. Measure points added before this will display a negative time in the summary.
    void reset(){start=std::chrono::high_resolution_clock::now();}
};

// /////////////////////////////

inline float fmodf_fast(float value,float divider)
{
    float ret;
    ret=value/divider;
    ret=ret-(long)ret;
    ret*=divider;
    return ret;
}

struct cosine_cache
{
    float values[91];

    cosine_cache()
    {
        for(int i=0;i<=90;i++)
            values[i]=cosf(i*3.14159265f/180.0f);
    }
};

struct cosine_cache_360
{
    float values[361];

    cosine_cache_360()
    {
        for(int i=0;i<=360;i++)
            values[i]=cosf(i*3.14159265f/180.0f);
    }
};

void SinCosf_cached_360(float degree, float* sin, float* cos)
{
    static cosine_cache_360 cache;
    if (degree < -360.0f)    // wrap input angle to -360...360
        degree = fmodf_fast(degree,360.0f);
    else if (degree >  360.0f)
        degree = fmodf_fast(degree,360.0f);

    float mod;
    float lower;
    float delta;
    int index;

    {   // sine
        float degree_sine=degree-90;
        if (degree_sine < 0.0f)      // translate for sine
            degree_sine += 360.0f;

        index=degree_sine;
        mod=degree_sine-index;
        lower=cache.values[index];
        delta=cache.values[index+1]-lower;
        *sin=lower+delta*mod;
    }
    {   // cosine
        if (degree < 0.0f)      // flip for cosine
            degree = -degree;

        index=degree;
        lower=cache.values[index];
        delta=cache.values[index+1]-lower;
        *cos=lower+delta*mod;
    }
}

// this is slow and currently nor working properly
void SinCosf_cached_90(float degree, float* sin, float* cos)
{
    static cosine_cache cache;
    if (degree < -360.0f)    // wrap input angle to -360...360
        degree = fmodf_fast(degree,360.0f);
    else if (degree >  360.0f)
        degree = fmodf_fast(degree,360.0f);

    float mod;
    float lower;
    float delta;
    int index;

    {   // sine
        float degree_sine=degree;
        if (degree_sine < 0.0f)      // translate for sine
            degree_sine += 180.0f;
        if (degree_sine > 90.0f)
        {
            degree_sine = 90.0f-degree;
            index=degree_sine;
            mod=degree_sine-index;
            lower=cache.values[90-index];
            delta=cache.values[89-index]-lower;
            *sin=lower+delta*mod;
        }
        else
        {
            index=degree_sine;
            mod=degree_sine-index;
            lower=cache.values[90-index];
            delta=cache.values[89-index]-lower;
            *sin=lower+delta*mod;
        }
    }
    {   // cosine
        if (degree < 0.0f)      // flip for cosine
            degree = -degree;
        if (degree > 90.0f)
        {
            degree = 90.0f-degree;
            lower=cache.values[index];
            delta=cache.values[index+1]-lower;
            *cos=-(lower+delta*mod);
        }
        else
        {
            lower=cache.values[index];
            delta=cache.values[index+1]-lower;
            *cos=lower+delta*mod;
        }
    }
}

void sincosf_fast3(float x, float* sin, float* cos)
{
    if (x < -6.28318531f)    // wrap input angle to -360..360
        x = fmodf_fast(x,6.28318531f);
    else if (x > 6.28318531f)
        x = fmodf_fast(x,6.28318531f);

    if (x < -3.14159265f)    // wrap input angle to -180..180
        x += 6.28318531f;
    else if (x > 3.14159265f)
        x -= 6.28318531f;

    const float pi = 3.14159265f;
    const float B = 4/pi;
    const float C = -4/(pi*pi);

    float y = B * x + C * x * fabsf(x);

    const float P = 0.225;

    *sin = P * (y * fabsf(y) - y) + y;   // Q * y + P * y * abs(y)

    x += 1.57079632f; // 90<=x<=270
    if (x > 3.14159265f)
        x -= 6.28318531f;
    y = B * x + C * x * fabsf(x);
    *cos = P * (y * fabsf(y) - y) + y;   // Q * y + P * y * abs(y)
}

void sincosf_fast2(float x, float* sin, float* cos)
{
    if (x < -6.28318531f)    // wrap input angle to -360..360
        x = fmodf_fast(x,6.28318531f);
    else if (x > 6.28318531f)
        x = fmodf_fast(x,6.28318531f);

    if (x < -3.14159265f)    // wrap input angle to -180..180
        x += 6.28318531f;
    else if (x > 3.14159265f)
        x -= 6.28318531f;

    float t;
    // compute sine and cosine: sin(x + 90°) = cos(x)
    if (x < 0)      // -180<=x<0
    {
        t = 1.27323954f * x + 0.405284735f * x * x;
        if (t < 0)
            *sin = 0.225f * (t *-t - t) + t;
        else
            *sin = 0.225f * (t * t - t) + t;

        x += 1.57079632f; // -90<=x<90
    }
    else            // 0<=x<=180
    {
        t = 1.27323954f * x - 0.405284735f * x * x;
        if (t < 0)
            *sin = 0.225f * (t *-t - t) + t;
        else
            *sin = 0.225f * (t * t - t) + t;

        x += 1.57079632f; // 90<=x<=270
        if (x > 3.14159265f)
            x -= 6.28318531f;
    }

    if (x < 0)
    {
        t = 1.27323954f * x + 0.405284735f * x * x;
        if (t < 0)
            *cos = 0.225f * (t *-t - t) + t;
        else
            *cos = 0.225f * (t * t - t) + t;
    }
    else
    {
        t = 1.27323954f * x - 0.405284735f * x * x;
        if (t < 0)
            *cos = 0.225f * (t *-t - t) + t;
        else
            *cos = 0.225f * (t * t - t) + t;
    }
}

void SinCosfFast(float x, float* sin, float* cos)
{
    x=fmodf_fast(x,360);

    // always wrap input angle to -180..180
    if (x < -180.0f)
        x += 360.0f;
    else
    if (x >  180.0f)
        x -= 360.0f;

    // compute sine
    if (x < 0)
        *sin = 0.02222222222f * x + 0.0001234567901f * x * x;
    else
        *sin = 0.02222222222f * x - 0.0001234567901f * x * x;

    //compute cosine: sin(x + PI/2) = cos(x)
    x += 90.0f;
    if (x >  180.0f)
        x -= 360.0f;

    if (x < 0)
        *cos = 0.02222222222f * x + 0.0001234567901f * x * x;
    else
        *cos = 0.02222222222f * x - 0.0001234567901f * x * x;
}

// based on the Low precision sine/cosine from http://lab.polygonal.de/?p=205. In my test it itself is ~4.5 times faster as the normal sinf.
void sincosf_fast(float x, float* sin, float* cos)
{
    x=fmodf_fast(x,3.14159265*2);

    // always wrap input angle to -180..180
    if (x < -3.14159265f)
        x += 6.28318531f;
    else
    if (x >  3.14159265f)
        x -= 6.28318531f;

    // compute sine
    if (x < 0)
        *sin = 1.27323954f * x + 0.405284735f * x * x;
    else
        *sin = 1.27323954f * x - 0.405284735f * x * x;

    //compute cosine: sin(x + PI/2) = cos(x)
    x += 1.57079632f;
    if (x >  3.14159265f)
        x -= 6.28318531f;

    if (x < 0)
        *cos = 1.27323954f * x + 0.405284735f * x * x;
    else
        *cos = 1.27323954f * x - 0.405284735f * x * x;
}

using namespace std;

void benchmark_fmodf()
{
    float f2;
    float f3=0;
    for(float f=-1000;f<1000;f+=55.5f)
    {
        cout<<f<<' ';
        f2=f/360.0f;
        f2=f2-(long)f2;
        f2*=360.0f;
        cout<<f2<<' '<<fmodf(f,360.0f)<<' '<<fmodf_fast(f,360.0f)<<endl;
    }

    for(int i=0;i<10;i++)
    {
        timer _("fmodf");
        for(int j=0;j<5;j++)
        for(float f=-2000;f<2000;f+=0.0001f)
            f3+=fmodf(f,360.0f);
    }
    for(int i=0;i<10;i++)
    {
        timer _("manual");
        for(int j=0;j<5;j++)
        for(float f=-2000;f<2000;f+=0.0001f)
        {
            f2=f/360.0f;
            f2=f2-(long)f2;
            f2*=360.0f;
            f3+=f2;
        }
    }
    for(int i=0;i<10;i++)
    {
        timer _("fmodf_fast");
        for(int j=0;j<5;j++)
        for(float f=-2000;f<2000;f+=0.0001f)
            f3+=fmodf_fast(f,360.0f);
    }
    cout<<f3<<endl;
}

int main()
{
    float sin;
    float cos;
    for(float i=-720;i<=720;i+=90.0f/4)
    {
        sincosf_fast2(i*3.14159265/180,&sin,&cos);
        std::cout<<"degree: "<<i<<" SinCosfFast2: \t"<<sin<<" \t"<<cos;
        sincosf_fast3(i*3.14159265/180,&sin,&cos);
        std::cout<<" \tSinCosfFast3 "<<sin<<" "<<cos;
        SinCosf_cached_90(i,&sin,&cos);
        std::cout<<" \tSinCosf_cached_90 "<<sin<<" "<<cos;
        SinCosf_cached_360(i,&sin,&cos);
        std::cout<<" \tSinCosf_cached_360 "<<sin<<" "<<cos;
        //SinCosfFast(i,&sin,&cos);
        //std::cout<<" F "<<" "<<sin<<" "<<cos;
        //sincosf_fast(i*3.14159265/180,&sin,&cos);
        //std::cout<<" _f "<<" "<<sin<<" "<<cos;
        sincosf(i*3.14159265/180,&sin,&cos);
        std::cout<<" \tsincosf: "<<sin<<" "<<cos<<std::endl;
        //std::cout<<i<<" "<<sinf(i*3.14159265/180)<<" "<<sinf_fast(i*3.14159265/180)<<" "<<sinf_fast2(i*M_PI/180)<<std::endl;
    }

    float s=0;
    float c=0;
    const float i_end=720.0f*2;
    {
        timer _("sinf cosf");
        for(float i=-720.0f;i<=i_end;i+=0.0001f)
        {
            s=sinf(i);
            c=cosf(i);
        }
    }
    {
        timer _("sincosf");
        for(float i=-720.0f;i<=i_end;i+=0.0001f)
            sincosf(i*3.14159265/180,&s,&c);
    }
    {
        timer _("sincosf_fast");
        for(float i=-720.0f;i<=i_end;i+=0.0001f)
            sincosf_fast(i*3.14159265/180,&s,&c);
    }
    {
        timer _("SinCosfFast");
        for(float i=-720.0f;i<=i_end;i+=0.0001f)
            SinCosfFast(i,&s,&c);
    }
    {
        timer _("sincosf_fast2 (precise and better branching)");
        for(float i=-720.0f;i<=i_end;i+=0.0001f)
            sincosf_fast2(i*3.14159265/180,&s,&c);
    }
    {
        timer _("sincosf_fast3 (precise and using fabsf instead of branches)");
        for(float i=-720.0f;i<=i_end;i+=0.0001f)
            sincosf_fast3(i*3.14159265/180,&s,&c);
    }
    {
        timer _("SinCosf_cached_90");
        for(float i=-720.0f;i<=i_end;i+=0.0001f)
            SinCosf_cached_90(i,&s,&c);
    }
    {
        timer _("SinCosf_cached_360");
        for(float i=-720.0f;i<=i_end;i+=0.0001f)
            SinCosf_cached_360(i,&s,&c);
    }
    std::cout<<s<<" "<<c<<std::endl;
}
