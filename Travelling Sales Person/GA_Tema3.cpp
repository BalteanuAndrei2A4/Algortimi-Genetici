#include <iostream>
#include <vector>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <random>
#include <chrono>
#include <cmath>
#include <string>
#include <cstring>
#include <fstream>
#include <map>
#include <math.h>

using namespace std;

#define POP_SIZE 100
#define PI 3.14159265358979323846  

int pop_size = 100;
int n = 403;
double minim_eval = 100000000;
vector <vector<int>> population;
vector <int> cromosome;
vector <int> unitate;
int distances[500][500]; int dimension = -1;
vector<double> fitnesses;
struct coordinates {
    double x;
    double y;
};
struct coordinates c[300];

void compute_distance(struct coordinates c[300], int nodes) // geographical distance
{
    for (int i = 0; i < nodes; i++)
        cout << i << " " << c[i].x << " " << c[i].y << endl;

    for (int i = 0; i < nodes - 1; i++)
    {
        for (int j = i + 1; j < nodes; j++)
        {
            double latitude_i = (PI * (floor(c[i].x) + 5 / 3 * (c[i].x - floor(c[i].x)))) / 180.0;
            double longitude_i = (PI * (floor(c[i].y) + 5 / 3 * (c[i].y - floor(c[i].y)))) / 180.0;
            //cout << "latitude of i " << latitude_i << " longitude of i " << longitude_i << endl;

            double latitude_j = (PI * (floor(c[j].x) + 5 / 3 * (c[j].x - floor(c[j].x)))) / 180.0;
            double longitude_j = (PI * (floor(c[j].y) + 5 / 3 * (c[j].y - floor(c[j].y)))) / 180.0;
            double R = 6378.388;

            double q1 = cos(longitude_i - longitude_j);
            double q2 = cos(latitude_i - latitude_j);
            double q3 = cos(latitude_i + latitude_j);
            // cout << q1 << " " << q2 << " " << q3 << endl;
            // cout<<acos( 0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3));

            distances[i][j] = floor(R * acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0);

            //cout << (value) << endl;
        }
    }
    /*
    for (int i = 0; i < nodes; i++)
    {
        for (int j = i + 1; j < nodes; j++)
        {
            cout<<distances[i][j]<<" ";
        }
        cout << endl;
    }
    */
}
void read_file_tsp(string filename)
{
    string line;
    ifstream infile;
    infile.open(filename);
    int count_line = 1;
    while (!infile.eof()) // to get all the lines.
    {

        getline(infile, line); // read && process the line 
        {
            int from_node = 0, x_from_node = 0, y_from_node = 0;
            int to_node, x_to_node, y_to_node;
            const char* temp = line.c_str();

            for (int i = 0; i < strlen(temp); i++)
            {

                if (temp[i] != ' ') // node
                {
                    int j = i;

                    if (from_node == 0)
                    {
                        while (temp[j] != ' ' && j < strlen(temp))
                        {
                            from_node *= 10;
                            from_node += temp[j] - '0';
                            j++;
                        }
                    }
                    else
                        if (x_from_node == 0)
                        {
                            int nat = 1;
                            float decimal = 0;
                            int cif = 0;
                            while (temp[j] != ' ' && j < strlen(temp))
                            {
                                if (temp[j] == '.')
                                {
                                    nat = 0;
                                    j++;
                                }
                                if (nat == 1)
                                {
                                    x_from_node *= 10;
                                    x_from_node += temp[j] - '0';
                                    j++;
                                }
                                else
                                {
                                    decimal *= 10;
                                    decimal += temp[j] - '0';
                                    cif++;
                                    j++;
                                }
                            }

                            if (decimal != 0)
                            {
                                c[from_node - 1].x = (x_from_node)+float(decimal) / pow(10, cif);
                            }
                            else
                            {
                                c[from_node - 1].x = x_from_node;
                            }

                        }
                        else
                            if (y_from_node == 0)
                            {
                                int nat = 1;
                                float decimal = 0;
                                int cif = 0;
                                while (temp[j] != ' ' && j < strlen(temp))
                                {
                                    //cout<<temp[j]<<endl;
                                    if (temp[j] == '.')
                                    {
                                        nat = 0;
                                        j++;
                                    }
                                    if (nat == 1)
                                    {
                                        y_from_node *= 10;
                                        y_from_node += temp[j] - '0';
                                        j++;
                                    }
                                    else
                                    {
                                        decimal *= 10;
                                        decimal += temp[j] - '0';
                                        cif++;
                                        j++;
                                    }

                                }
                                if (decimal != 0)
                                {
                                    c[from_node - 1].y = y_from_node + decimal / pow(10, cif);
                                }
                                else
                                {
                                    c[from_node - 1].y = y_from_node;
                                }
                            }

                    i = j;
                }
            }
            count_line++;
        }
    }
    infile.close();
    count_line--;
    compute_distance(c, count_line);
    /*
    for (int i = 0; i < count_line; i++)
        cout << i << " " << c[i].x << " " << c[i].y << endl;
        */
        /*
        for (int i = 0; i < count_line - 1; i++)
        {
            for (int j = i + 1; j < count_line; j++)
            {
                distances[i][j] = int(sqrt((c[i].x - c[j].x) * (c[i].x - c[j].x) + (c[i].y - c[j].y) * (c[i].y - c[j].y)));
                //cout << distances[i][j] << " ";
            }
           // cout << endl;
        }
        */
}

int lineno = -1;
void read_file_atsp(string filename)
{
    string line;
    ifstream infile;
    infile.open(filename);
    int count_line = 1;
    int count_nodes = 0;
    while (!infile.eof()) // to get all the lines.
    {

        getline(infile, line); // read && process the line 

        const char* temp = line.c_str(); // cout << line << endl;
        if (dimension == -1) // dimension of matrix
        {
            lineno++;
            dimension = atoi(line.c_str());
            // cout << dimension << endl;
        }
        else
        {
            for (int ii = 0; ii < strlen(temp); ii++)
            {

                int weight = 0;
                if (temp[ii] != ' ') // node
                {
                    while (temp[ii] != ' ' && ii < strlen(temp))
                    {
                        weight *= 10;
                        weight += temp[ii] - '0';
                        ii++;
                    }

                    // cout << "weight from file "<<weight << endl;
                    // cout << "distances[" << count_nodes << "][" << lineno << "] = " << weight << endl;
                    distances[lineno][count_nodes] = weight;
                    count_nodes++;

                    if (count_nodes == dimension)
                    {
                        count_nodes = 0;
                        lineno++;
                    }
                }
            }

        }
    }
}

vector<vector<int>> generate_population() {
    vector<vector<int>> result;
    int max_random;
    for (int i = 0; i < pop_size; ++i) {
        max_random = dimension;
        vector<int> row(dimension - 1);
        for (int& x : row) {
            x = rand() % max_random;
            max_random--;
        }
        result.push_back(row);
        /*
        for (auto ii = 0; ii < row.size(); ii++)
        {
            cout << row.at(ii) << " ";
        }
        cout << endl;
        */
    }

    return result;
}

double random_double(double range_from, double range_to) {
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_real_distribution<> distr(range_from, range_to);
    return distr(generator);
}

int random_int(int range_from, int range_to) {
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_int_distribution<> distr(range_from, range_to);
    return distr(generator);
}

double evaluate(vector <int> vec) {

    double result = 0;
    for (int i = 0; i < vec.size() - 1; ++i)
    {
        result += distances[vec[i]][vec[i + 1]];

    }

    result += distances[vec[0]][vec[vec.size() - 1]];

    //cout << "from evaluate function " << int(result) << endl;
    return result;
}

double fitness_function(vector<int> vec) {
    double result;
    result = evaluate(vec);
    double result_fitness;
    result_fitness = 1 / (result);

    return result_fitness;
}

void mutate(double pm, vector<vector<int>>& v) {
    double r;
    int bit;
    int max_random;
    for (int i = 0; i < pop_size; ++i) {
        max_random = n;
        for (int j = 0; j < n - 1; j++) {
            r = random_double(0.0, 1.0);
            if (r < pm) {
                v[i][j] = (v[i][j] + 1) % max_random;
            }
            max_random--;
        }
    }
}

vector <int> unitate_construct()
{
    vector <int> result(dimension);
    for (int i = 0; i < dimension; ++i)
    {
        result[i] = i;
    }
    return result;
}

vector <int> delete_poz(int poz, vector<int> v)
{
    int n = v.size();
    int j = 0;
    vector <int> result(n - 1);
    for (int i = 0; i < result.size() + 1; ++i) {

        if (poz != i)
            result[j++] = v[i];
    }
    return result;
}

vector <int> decode(vector <int> cromosome)
{
    vector <int> result(dimension);
    unitate = unitate_construct();
    /*
    cout << "Unitate contine : \n";
    for (int i = 0; i < dimension; ++i)
    {
        cout<<unitate[i]<<" ";
    }
    cout << "\n";
   */ // cout << "done\n";

    int i;
    for (i = 0; i < cromosome.size(); ++i) {
        result[i] = unitate[cromosome[i]];
        unitate = delete_poz(cromosome[i], unitate);
    }
    result[i] = unitate[0];
    /*for (int i = 0; i < dimension; ++i)
    {
        cout << result[i] << " ";
    }
    cout << "\n";
    */
    return result;
}

vector <vector<int>> selection(vector <vector<int>> v)
{
    //cout << "selection called\n";
    vector <double> evaluate_vector;
    vector <double> individual_selection;
    vector <vector<int>> result;

    for (int i = 0; i < population.size(); i++)
    {
        population[i] = decode(population[i]);
    }

    int k = 0;
    double r, value;
    double Total_fitness = 0;

    for (int j = 0; j < pop_size; ++j) {

        value = fitness_function(decode(v[j]));
        evaluate_vector.push_back(value);
        Total_fitness += evaluate_vector[j];
    }

    for (int j = 0; j < pop_size; ++j) {
        individual_selection.push_back(evaluate_vector[j] / Total_fitness);
    }

    double accumulated_selection[POP_SIZE + 1];
    accumulated_selection[0] = 0;

    for (int j = 0; j < pop_size; ++j) {
        accumulated_selection[j + 1] = accumulated_selection[j] + individual_selection[j];
    }

    int p, ok;

    for (int i = 0; i < pop_size; ++i)
    {
        r = random_double(0.0, 1.0);
        ok = 1;
        for (int j = 0; j < pop_size; ++j) {
            if (accumulated_selection[j] <= r && r <= accumulated_selection[j + 1]) {
                result.push_back(v[j]);
                break;
            }
        }
    }
    
  
     
    fitnesses = evaluate_vector;
    return result;
}

void cross_over_switch(int x, int y, vector<vector<int>>& v) {
    vector <int> c1;
    vector <int> c2;

    int poz = random_int(0, dimension - 2);
    // cout << poz << endl;
    for (int j = 0; j < poz; ++j)
        c1.push_back(v[x][j]);
    for (int j = poz; j < v[x].size(); ++j)
        c1.push_back(v[y][j]);

    for (int j = 0; j < poz; ++j)
        c2.push_back(v[y][j]);
    for (int j = poz; j < v[x].size(); ++j)
        c2.push_back(v[x][j]);

    for (int i = 0; i < c1.size(); i++)
    {
        v[x][i] = c1[i];
    }
    for (int i = 0; i < c2.size(); i++)
    {
        v[y][i] = c2[i];
    }
}

void cross_over(double pcx, vector <vector<int>>& v) {
    struct
    {
        vector <double> probability;
        int* indice = new int[pop_size];
    }lista;

    double r;

    for (int i = 0; i < pop_size; ++i) {
        r = random_double(0.0, 1.0);
        lista.probability.push_back(r);
        lista.indice[i] = i;
    }

    for (int i = 0; i < pop_size - 1; ++i) {
        for (int j = i + 1; j < pop_size; ++j) {
            if (lista.probability[i] > lista.probability[j])
            {
                double aux;
                aux = lista.probability[i];
                lista.probability[i] = lista.probability[j];
                lista.probability[j] = aux;

                int aux_indice;
                aux_indice = lista.indice[i];
                lista.indice[i] = lista.indice[j];
                lista.indice[j] = aux_indice;
            }
        }
    }

    for (int i = 0; i < pop_size - 2; i = i + 2) {
        if (lista.probability[i] <= pcx && lista.probability[i + 1] <= pcx) {
            // cout << "Facem un cross-over " << lista.indice[i] << " cu " << lista.indice[i + 1] << endl;
            cross_over_switch(lista.indice[i], lista.indice[i + 1], v);
        }
    }
}

double GA()
{
    population = generate_population();
    int optim = evaluate(decode((population.at(0))));
    int generation = 0;
    vector <vector<int>> next_gen;
    vector<int> aux;
    while (generation < 1000)
    {
        generation++;
        population = selection(population);
        mutate(0.01, population);
        cross_over(0.2, population);
       
    }
 for (int i = 0; i < pop_size; i++)
        {
            int temp_sol = evaluate(decode(population.at(i)));
            if (temp_sol < optim)
            {
                optim = temp_sol;
                aux = decode(population.at(i));
            }
        }
     /* 
     cout << "PERMUTAREA\n";
     sort(aux.begin(), aux.end());
     for (int i = 0; i < aux.size(); i++)
     {
         cout << aux[i] << " ";
     }
   */
    return optim;
}

void improve(vector<int> candidate)
{

    double r;
    int bit;
    int max_random;
    int pos = n;

    for (int j = 0; j < n - 1; j++)
    {
        candidate[j] = (candidate[j] + 1) % pos;
        pos--;
    }
}

vector <int> neighbor(int poz, vector <int> v) {
    int max_random = n;
    vector <int> result(n - 1);
    for (int i = 0; i < v.size(); ++i)
    {
        if (i == poz)
            result[i] = (v[i] + 1) % max_random;
        else result[i] = v[i];
        max_random--;
    }
    return result;
}
vector <int> generate_random() {
    vector <int> result;
    int max_random = n;
    for (int i = 0; i < n - 1; i++) {
        int x = rand() % max_random;
        result.push_back(x);
        max_random--;
    }
    return result;
}

double SA()
{
    int t = 0;
    int T = 100;
    vector <int> vc;
    vector <int> vn;
    double r;
    double best = 100000000;
    vc = generate_random();
    double eval_vc;
    eval_vc = evaluate(decode(vc));
    do {
        int vec = 0;
        int k = 0;
        do {
            k++;
            r = ((double)rand() / (RAND_MAX));
            vn = neighbor(vec, vc);
            vec++;
            double best_n = evaluate(decode(vn));
            vector <int> best_array_n;
            best_array_n = vn;

            while (vec < vc.size())
            {
                vn = neighbor(vec, vc);
                vec++;
                if (evaluate(decode(vn)) < evaluate(decode(best_array_n)))
                {
                    best_array_n = vn;
                }
            }
            if (evaluate(decode(best_array_n)) < evaluate(decode(vc)))
            {
                vc = best_array_n;
                vec = 0;
            }
            else if (r < exp(-abs(evaluate(decode(best_array_n)) - evaluate(decode(vc))) / T))
            {
                vc = best_array_n;
                vec = 0;
            }
        } while (k < 1000);

        T = T * 0.99;
        if (evaluate(decode(vc)) < best)
            best = evaluate(decode(vc));

    } while (T > pow(10, -8));
    return best;

}


int main()
{

   //read_file_atsp("br17_atsp.txt");
   
   cout << " GENETIC ALGORTIHM \n";
   
    for (int i = 0; i < 30; i++)
    {
      cout << GA()<<endl;
    }
    cout << endl;
   
    
    cout << " SIMMULATED ANNEALING \n";
    for (int i = 0; i < 30; i++)
    {
        cout  << SA() << endl;
        
    }
    
    return 0;
}