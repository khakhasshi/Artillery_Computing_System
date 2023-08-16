#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_AMMO 10
#define FILE_NAME "ammo.txt"
typedef struct {
    char name[20];
    double initialVelocity;
    double mass;
    double dragCoefficient;
    double crossSectionalArea;
    double explosionRadius;
    double maxRange; // 平地最大射程
} Ammo;

Ammo ammoList[MAX_AMMO];
int ammoCount = 0;

void saveAmmoToFile();
void loadAmmoFromFile();
double calculateDistance(double angle, double initialVelocity, double dragCoefficient, double crossSectionalArea, double mass);
double calculateAirResistance(double velocity, double dragCoefficient, double crossSectionalArea);
double calculateMaxRange(double initialVelocity, double dragCoefficient, double crossSectionalArea, double mass);
double bisectIterate(double targetX, double targetY, double initialVelocity, double dragCoefficient, double crossSectionalArea, double mass);
void shoot();
void addAmmo();
void viewAmmo();
void deleteAmmo(); 


int main() {
    loadAmmoFromFile();

    FILE *file = fopen(FILE_NAME, "r");
    if (file == NULL) {
        file = fopen(FILE_NAME, "w");
        if (file != NULL) {
            fclose(file);
            printf("Created %s file.\n", FILE_NAME);
        } else {
            printf("Error creating %s file.\n", FILE_NAME);
        }
    } else {
        fclose(file);
    }

     int choice;
    while (1) {
        printf("\nFire Cannon Management System\n");
        printf("1. Shoot\n");
        printf("2. Add Ammo\n");
        printf("3. Delete Ammo\n"); // New option for deleting ammo
        printf("4. View Ammo\n");
        printf("5. Exit\n");
        printf("Enter your choice: ");
        scanf("%d", &choice);

        switch (choice) {
            case 1:
                shoot();
                break;
            case 2:
                addAmmo();
                break;
            case 3:
                deleteAmmo(); // Call the new function for deleting ammo
                break;
            case 4:
                viewAmmo();
                break;
            case 5:
                saveAmmoToFile();
                printf("Exiting program.\n");
                exit(0);
            default:
                printf("Invalid choice. Please try again.\n");
                break;
        }
    }

    return 0;
}

// ...（前面的代码不变）

void saveAmmoToFile() {
    FILE *file = fopen(FILE_NAME, "w");
    if (file == NULL) {
        printf("Error opening file.\n");
        return;
    }

    fprintf(file, "%d\n", ammoCount);
    for (int i = 0; i < ammoCount; i++) {
        fprintf(file, "%s %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf\n",
                ammoList[i].name, ammoList[i].initialVelocity, ammoList[i].mass,
                ammoList[i].dragCoefficient, ammoList[i].crossSectionalArea,
                ammoList[i].explosionRadius, ammoList[i].maxRange);
    }

    fclose(file);
}

void deleteAmmo() {
    if (ammoCount == 0) {
        printf("No ammo available to delete.\n");
        return;
    }

    printf("Enter the index of the ammo to delete: ");
    int indexToDelete;
    scanf("%d", &indexToDelete);

    if (indexToDelete < 1 || indexToDelete > ammoCount) {
        printf("Invalid index.\n");
        return;
    }

    indexToDelete--; // Adjust for 0-based index

    for (int i = indexToDelete; i < ammoCount - 1; i++) {
        ammoList[i] = ammoList[i + 1];
    }

    ammoCount--;

    // Save the updated ammo list to the file
    saveAmmoToFile();

    printf("Ammo deleted successfully.\n");
}


void loadAmmoFromFile() {
    FILE *file = fopen(FILE_NAME, "r");
    if (file == NULL) {
        printf("No ammo file found.\n");
        return;
    }

    fscanf(file, "%d", &ammoCount);
    for (int i = 0; i < ammoCount; i++) {
        fscanf(file, "%s %lf %lf %lf %lf %lf %lf %lf",
            ammoList[i].name, &ammoList[i].initialVelocity, &ammoList[i].mass,
            &ammoList[i].dragCoefficient, &ammoList[i].crossSectionalArea,
            &ammoList[i].explosionRadius, &ammoList[i].maxRange);
    }

    fclose(file);
}

// ...（后面的代码不变）


double calculateDistance(double angle, double initialVelocity, double dragCoefficient, double crossSectionalArea, double mass) {
    double time = 0.0;
    double dt = 0.001;
    double velocityX = initialVelocity * cos(angle * M_PI / 180);
    double velocityY = initialVelocity * sin(angle * M_PI / 180);
    double positionX = 0.0;
    double positionY = 0.0;
    double distance = 0.0;
    double accelerationY = -9.81;

    while (positionY >= 0) {
        double airResistance = calculateAirResistance(sqrt(pow(velocityX, 2) + pow(velocityY, 2)), dragCoefficient, crossSectionalArea);
        double accelerationX = -airResistance * velocityX / sqrt(pow(velocityX, 2) + pow(velocityY, 2));
        velocityX += accelerationX * dt;
        velocityY += accelerationY * dt;
        positionX += velocityX * dt;
        positionY += velocityY * dt;
        distance += sqrt(pow(velocityX, 2) + pow(velocityY, 2)) * dt;
        time += dt;
    }

    return distance;
}

double calculateAirResistance(double velocity, double dragCoefficient, double crossSectionalArea) {
    double airDensity = 1.225;

    return 0.5 * dragCoefficient * airDensity * pow(velocity, 2) * crossSectionalArea;
}

double calculateMaxRange(double initialVelocity, double dragCoefficient, double crossSectionalArea, double mass) {
    double maxRange = 0.0;
    double increment = 1.0;
    double angle = 45.0;
    double distance;

    while (angle >= 0.0) {
        distance = calculateDistance(angle, initialVelocity, dragCoefficient, crossSectionalArea, mass);
        if (distance > maxRange) {
            maxRange = distance;
        }
        angle -= increment;
    }

    return maxRange;
}

double bisectIterate(double targetX, double targetY, double initialVelocity, double dragCoefficient, double crossSectionalArea, double mass) {
    double lowerAngle = 0.0;
    double upperAngle = 45.0;
    double error = 1.0; // Initial error value
    double tolerance = 0.001; // Convergence tolerance
    int maxIterations = 100; // Maximum number of iterations

    int iteration = 0;
    double angle = 0.0;

    while (fabs(error) > tolerance && iteration < maxIterations) {
        angle = (lowerAngle + upperAngle) / 2.0;
        double distance = calculateDistance(angle, initialVelocity, dragCoefficient, crossSectionalArea, mass);
        error = distance - sqrt(pow(targetX, 2) + pow(targetY, 2));

        if (error < 0) {
            lowerAngle = angle;
        } else {
            upperAngle = angle;
        }

        iteration++;
    }

    if (angle <= 0 || angle >= 90) {
        printf("Unable to estimate angle. Target is out of range.\n");
        return -1.0;
    }

    return angle;
}

void shoot() {
    int ammoIndex;
    double targetX, targetY;

    printf("\nAvailable Ammo:\n");
    for (int i = 0; i < ammoCount; i++) {
        printf("%d. %s (Max Range: %.2lf)\n", i + 1, ammoList[i].name, ammoList[i].maxRange);
    }

    printf("Enter the ammo index: ");
    scanf("%d", &ammoIndex);
    ammoIndex--;

    if (ammoIndex < 0 || ammoIndex >= ammoCount) {
        printf("Invalid ammo index.\n");
        return;
    }

    printf("Enter the target coordinates (x, y): ");
    scanf("%lf %lf", &targetX, &targetY);

    double maxRange = ammoList[ammoIndex].maxRange;
    double distanceToTarget = sqrt(pow(targetX, 2) + pow(targetY, 2));

    if (distanceToTarget > maxRange) {
        printf("Target is out of range.\n");
        return;
    }

    double initialVelocity = ammoList[ammoIndex].initialVelocity;
    double dragCoefficient = ammoList[ammoIndex].dragCoefficient;
    double crossSectionalArea = ammoList[ammoIndex].crossSectionalArea;
    double mass = ammoList[ammoIndex].mass;

    double angle = bisectIterate(targetX, targetY, initialVelocity, dragCoefficient, crossSectionalArea, mass);

    if (angle != -1.0 && angle >= 0.0 && angle <= 45.0) {
        printf("Estimated shooting angle: %.12lf degrees\n", angle);
        double mils = angle * 17.777777777778;
        printf("NATO MILS: %.12lf mils\n", mils);

        // Simulate projectile motion with air resistance using Runge-Kutta method
        double velocity = initialVelocity;
        double velocityX = velocity * cos(angle * M_PI / 180);
        double velocityY = velocity * sin(angle * M_PI / 180);
        double dt = 0.00001;
        double positionX = 0.0;
        double positionY = 0.0;
        double distance = 0.0;
        double accelerationY = -9.81;
        double airResistance;

        while (positionY >= 0) {
            airResistance = calculateAirResistance(velocity, dragCoefficient, crossSectionalArea);

            // Runge-Kutta method - Step 1
            double accelerationX1 = -airResistance * velocityX / velocity;
            double accelerationY1 = accelerationY;

            // Runge-Kutta method - Step 2
            double velocityX2 = velocityX + accelerationX1 * dt / 2.0;
            double velocityY2 = velocityY + accelerationY1 * dt / 2.0;
            double accelerationX2 = -calculateAirResistance(sqrt(pow(velocityX2, 2) + pow(velocityY2, 2)), dragCoefficient, crossSectionalArea) * velocityX2 / velocity;
            double accelerationY2 = accelerationY;

            // Runge-Kutta method - Step 3
            double velocityX3 = velocityX + accelerationX2 * dt / 2.0;
            double velocityY3 = velocityY + accelerationY2 * dt / 2.0;
            double accelerationX3 = -calculateAirResistance(sqrt(pow(velocityX3, 2) + pow(velocityY3, 2)), dragCoefficient, crossSectionalArea) * velocityX3 / velocity;
            double accelerationY3 = accelerationY;

            // Runge-Kutta method - Step 4
            double velocityX4 = velocityX + accelerationX3 * dt;
            double velocityY4 = velocityY + accelerationY3 * dt;
            double accelerationX4 = -calculateAirResistance(sqrt(pow(velocityX4, 2) + pow(velocityY4, 2)), dragCoefficient, crossSectionalArea) * velocityX4 / velocity;
            double accelerationY4 = accelerationY;

            // Runge-Kutta method - Update variables
            velocityX += (accelerationX1 + 2.0 * accelerationX2 + 2.0 * accelerationX3 + accelerationX4) * dt / 6.0;
            velocityY += (accelerationY1 + 2.0 * accelerationY2 + 2.0 * accelerationY3 + accelerationY4) * dt / 6.0;
            velocity = sqrt(pow(velocityX, 2) + pow(velocityY, 2));
            positionX += velocityX * dt;
            positionY += velocityY * dt;
            distance += velocity * dt;
        }

        // Calculate kinetic energy at target
        double kineticEnergy = 0.5 * mass * pow(velocity, 2);

        printf("Velocity at target: %.12lf m/s\n", velocity);
        printf("Kinetic energy at target: %.12lf J\n", kineticEnergy);
    } else {
        printf("Unable to estimate angle. Target is out of range.\n");
    }
}


void addAmmo() {
    if (ammoCount == MAX_AMMO) {
        printf("Cannot add more ammo. Maximum limit reached.\n");
        return;
    }

    printf("\nEnter the ammo details:\n");
    printf("Name: ");
    scanf("%s", ammoList[ammoCount].name);
    printf("Initial Velocity: ");
    scanf("%lf", &ammoList[ammoCount].initialVelocity);
    printf("Mass: ");
    scanf("%lf", &ammoList[ammoCount].mass);
    printf("Drag Coefficient: ");
    scanf("%lf", &ammoList[ammoCount].dragCoefficient);
    printf("Cross-sectional Area: ");
    scanf("%lf", &ammoList[ammoCount].crossSectionalArea);
    printf("Explosion Radius: ");
    scanf("%lf", &ammoList[ammoCount].explosionRadius);

    // Calculate the max range for the added ammo
    ammoList[ammoCount].maxRange = calculateMaxRange(ammoList[ammoCount].initialVelocity,
                                                     ammoList[ammoCount].dragCoefficient,
                                                     ammoList[ammoCount].crossSectionalArea,
                                                     ammoList[ammoCount].mass);

    ammoCount++;

    // 保存弹药信息到文件
    FILE *file = fopen(FILE_NAME, "a");
    if (file == NULL) {
        printf("Error opening file for appending.\n");
        return;
    }

    fprintf(file, "%s %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf\n",
            ammoList[ammoCount - 1].name, ammoList[ammoCount - 1].initialVelocity,
            ammoList[ammoCount - 1].mass, ammoList[ammoCount - 1].dragCoefficient,
            ammoList[ammoCount - 1].crossSectionalArea, ammoList[ammoCount - 1].explosionRadius,
            ammoList[ammoCount - 1].maxRange);

    fclose(file);

    printf("Ammo added successfully.\n");
}
void viewAmmo() {
    if (ammoCount == 0) {
        printf("No ammo available.\n");
        return;
    }

    printf("\nAvailable Ammo:\n");
    for (int i = 0; i < ammoCount; i++) {
        printf("Name: %s\n", ammoList[i].name);
        printf("Initial Velocity: %.2lf\n", ammoList[i].initialVelocity);
        printf("Mass: %.2lf\n", ammoList[i].mass);
        printf("Drag Coefficient: %.2lf\n", ammoList[i].dragCoefficient);
        printf("Cross-sectional Area: %.2lf\n", ammoList[i].crossSectionalArea);
        printf("Explosion Radius: %.2lf\n", ammoList[i].explosionRadius);
        printf("Max Range: %.2lf\n", ammoList[i].maxRange);
        printf("\n");
    }
}

